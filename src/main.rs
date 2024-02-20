// const CLOCK_RATE: HertzU32 = 1323000;
const PDM_DECIMATION: u8 = 64;
const SAMPLE_RATE: u32 = 32000 * 64;
const PI: f32 = 3.14159;
const SINCN: u8 = 3;
const FILTER_GAIN: u8 = 16;
const MAX_VOLUME: u8 = 64;
mod pdm;
use pdm::PdmFilters;

pub mod cic;
pub mod fir;
struct PDMFilter {
    lut: Vec<u32>,
    fs: u32,
    coef: [u32; 3],
    sub_const: u32,
    old_out: i64,
    old_in: i64,
    oldz: i64,
    hp_alpha: u32,
    lp_alpha: u32,
    div_const: u32,
}
impl PDMFilter {
    pub fn new() -> PDMFilter {
        PDMFilter {
            lut: Vec::with_capacity((SINCN * PDM_DECIMATION / 8) as usize * 256),
            fs: SAMPLE_RATE,
            coef: [0, 0, 0],
            sub_const: 0,
            old_out: 0,
            old_in: 0,
            oldz: 0,
            lp_alpha: 0,
            hp_alpha: 0,
            div_const: 0,
        }
    }
    pub fn init(&mut self) {
        let lp_hz: f32 = SAMPLE_RATE as f32 / 128.0;
        let hp_hz: f32 = 10.0;
        if lp_hz != 0.0 {
            self.lp_alpha = (lp_hz as f32 * 256.0 / (lp_hz + self.fs as f32 / (2.0 * PI))) as u32;
        }
        if hp_hz != 0.0 {
            self.hp_alpha = (self.fs as f32 * 256.0 / (2.0 * PI * hp_hz + self.fs as f32)) as u32;
        }
        println!(
            "LP {} alpha {} HP {} alpha {}",
            lp_hz, hp_hz, self.lp_alpha, self.hp_alpha
        );
        let mut sinc = Vec::with_capacity(PDM_DECIMATION as usize);
        for _ in 0..PDM_DECIMATION {
            sinc.push(1)
        }
        let sinc2 = convolve(&sinc, &sinc);
        let sinc2 = convolve(&sinc2, &sinc);
        let mut sum: u32 = 0;
        let mut coef: Vec<u32> = Vec::with_capacity((PDM_DECIMATION * SINCN + 2) as usize);
        coef.push(0);
        for j in 0..SINCN {
            for i in 0..PDM_DECIMATION {
                // seems like this shouldn't happen...
                if (j * PDM_DECIMATION + i) as usize >= sinc2.len() {
                    coef.push(0);
                } else {
                    coef.push(sinc2[(j * PDM_DECIMATION + i) as usize] as u32);
                    sum += sinc2[(j * PDM_DECIMATION + i) as usize] as u32;
                }
                // if j >= 1 {
                // break;?
                // }
            }
        }
        coef.push(0);
        coef.push(0);
        self.sub_const = sum >> 1;
        self.div_const = self.sub_const * MAX_VOLUME as u32 / 32768 / FILTER_GAIN as u32;
        if self.div_const == 0 {
            self.div_const = 1;
        }

        println!(
            "Sum {} DIv const {} sub const {}",
            sum, self.div_const, self.sub_const
        );
        for s in 0..SINCN {
            let offset: usize = (s * PDM_DECIMATION) as usize;

            for c in 0..256 {
                for d in 0..(PDM_DECIMATION / 8) as usize {
                    self.lut.push(
                        (c >> 7) * coef[offset + d * 8]
                            + ((c >> 6) & 0x01) * coef[offset + d * 8 + 1]
                            + ((c >> 5) & 0x01) * coef[offset + d * 8 + 2]
                            + ((c >> 4) & 0x01) * coef[offset + d * 8 + 3]
                            + ((c >> 3) & 0x01) * coef[offset + d * 8 + 4]
                            + ((c >> 2) & 0x01) * coef[offset + d * 8 + 5]
                            + ((c >> 1) & 0x01) * coef[offset + d * 8 + 6]
                            + ((c) & 0x01) * coef[offset + d * 8 + 7],
                    );
                    // if s == 0 && c < 5 {
                    //     println!(
                    //         "s={} c={} d={} index {} =  {}",
                    //         s,
                    //         c,
                    //         d,
                    //         self.lut.len() - 1,
                    //         self.lut[self.lut.len() - 1]
                    //     );
                    // }
                }
            }
        }
        println!("lut is {}", self.lut.len());
    }

    fn filter(&mut self, data: &[u8], volume: u32) -> Vec<u16> {
        let mut dataout: Vec<u16> = Vec::with_capacity(data.len());

        let mut old_out: i64 = self.old_out;
        let mut old_in: i64 = self.old_in;
        let mut oldz: i64 = self.oldz;
        println!("Filtering {}", data.len());
        // skip PDM_DECIMATION bits, so change to bytes to match data u8 type
        for i in (0..data.len() - 8).step_by(PDM_DECIMATION as usize >> 3) {
            // println!("Starting {} stopping at {}", i, data.len() - 8);
            let index = i as usize;

            // 3 stages?
            let z0 = filter_table_mono_64(&self.lut, &data[index..index + 8], 0);
            let z1 = filter_table_mono_64(&self.lut, &data[index..index + 8], 1);
            let z2 = filter_table_mono_64(&self.lut, &data[index..index + 8], 2);
            // println!("z0 {} z1 {} z2 {}", z0, z1, z2);

            let mut z: i64 = self.coef[1] as i64 + z2 as i64 - self.sub_const as i64;
            // println!("z now is {}", z);
            self.coef[1] = self.coef[0] + z1 as u32;
            self.coef[0] = z0 as u32;

            old_out = (self.hp_alpha as i64 * (old_out + z - old_in)) >> 8;
            old_in = z;

            oldz = ((256 - self.lp_alpha as i64) * oldz + self.lp_alpha as i64 * old_out) >> 8;
            z = oldz * volume as i64;
            // println!("z now is {}", z);

            z = round_div(z, self.div_const as i64);

            z = satural_lh(z, -32700 as i64, 32700 as i64);
            // println!("Z is {}", z);
            // how is this right if satural lh can return -32700
            dataout.push(z as u16);
        }
        self.old_out = old_out;
        self.old_in = old_in;
        self.oldz = oldz;
        return dataout;
    }
}

// round(a/b)
fn round_div(a: i64, b: i64) -> i64 {
    if a > 0 {
        return (a + b / 2) / b;
    }
    return (a - b / 2) / b;
}
// clip???
fn satural_lh(n: i64, l: i64, h: i64) -> i64 {
    if n < l {
        return l;
    } else if n > h {
        return h;
    }
    return n;
}

// apply weights on each bit of input data
fn filter_table_mono_64(lut: &Vec<u32>, data: &[u8], s: u8) -> u32 {
    // let val = lut[(s * PDM_DECIMATION + data[0] * PDM_DECIMATION / 8) as usize];
    // println!(
    //     "Getting data at {} index {}",
    //     data[0],
    //     (data[0] as usize * PDM_DECIMATION as usize / 8)
    // );
    // println!(
    //     "Getting data at {} index {}",
    //     data[1],
    //     (data[1] as usize * PDM_DECIMATION as usize / 8 + 1)
    // );
    // println!(
    //     "Getting data at {} index {}",
    //     data[2],
    //     (data[2] as usize * PDM_DECIMATION as usize / 8 + 2)
    // );
    // println!(
    //     "Getting data at {} index {}",
    //     data[3],
    //     (data[3] as usize * PDM_DECIMATION as usize / 8 + 3)
    // );
    return lut[(data[0] as usize * PDM_DECIMATION as usize / 8) as usize]
        + lut[(data[1] as usize * PDM_DECIMATION as usize / 8 + 1) as usize]
        + lut[(data[2] as usize * PDM_DECIMATION as usize / 8 + 2) as usize]
        + lut[(data[3] as usize * PDM_DECIMATION as usize / 8 + 3) as usize]
        + lut[(data[4] as usize * PDM_DECIMATION as usize / 8 + 4) as usize]
        + lut[(data[5] as usize * PDM_DECIMATION as usize / 8 + 5) as usize]
        + lut[(data[6] as usize * PDM_DECIMATION as usize / 8 + 6) as usize]
        + lut[(data[7] as usize * PDM_DECIMATION as usize / 8 + 7) as usize];
}
fn convolve(signal: &Vec<u16>, kernel: &Vec<u16>) -> Vec<u16> {
    let outlen = signal.len() + kernel.len() - 1;
    let mut out: Vec<u16> = Vec::with_capacity(outlen);

    for n in 0..outlen {
        let mut kmin: usize = 0;
        if n >= kernel.len() - 1 {
            kmin = n - (kernel.len() - 1)
        }

        let mut kmax: usize = signal.len() - 1;
        if n < signal.len() - 1 {
            kmax = n
        }
        let mut acc: u16 = 0;

        for k in kmin..=kmax {
            acc += signal[k] * kernel[n - k]
        }

        out.push(acc);
    }
    return out;
}
use std::fs;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::{mem, slice};

fn get_file_as_byte_vec(filename: &String) -> Vec<u8> {
    let mut f = File::open(&filename).expect("no file found");
    let metadata = fs::metadata(&filename).expect("unable to read metadata");
    let mut buffer = vec![0; metadata.len() as usize];
    f.read(&mut buffer).expect("buffer overflow");

    buffer
}
fn main() {
    let mut filter = PDMFilter::new();
    filter.init();
    let filename = String::from("test_audio.raw");
    let buffer = get_file_as_byte_vec(&filename);

    // for val in &buffer[0..200] {
    //     print!("Val is {}", val);
    // }
    // println!("Have {}", buffer.len());
    // let filtered = filter.filter(&buffer[0..], 64);
    // println!("Filtered data {}", filtered.len());
    // for val in &filtered[..10] {
    //     println!("PCM Val is {}", val);
    // }
    // let slice_u8: &[u8] = unsafe {
    //     slice::from_raw_parts(
    //         filtered.as_ptr() as *const u8,
    //         filtered.len() * mem::size_of::<u16>(),
    //     )
    // };
    // println!("Filtered u8 data {}", slice_u8.len());

    // match fs::write("test_audio.pcm", &slice_u8) {
    //     Ok(()) => {
    //         println!(
    //             "Saved Audio file {} bytes are {}",
    //             "test_audio.pcm",
    //             slice_u8.len()
    //         );
    //     }
    //     Err(e) => {
    //         println!(
    //             "Failed writing Audio file to storage at {}, reason: {}",
    //             "test_audio.cpm", e
    //         );
    //     }
    // }

    let mut out: [f32; 80000] = [0f32; 80000];
    let mut pdm_processor = PdmFilters::new();
    pdm_processor.process_pdm(&buffer[..], 1, &mut out[..]);
    let slice_u8: &[u8] = unsafe {
        slice::from_raw_parts(out.as_ptr() as *const u8, out.len() * mem::size_of::<f32>())
    };
    match fs::write("test_audio.pcm", &slice_u8) {
        Ok(()) => {
            println!(
                "Saved Audio file {} bytes are {}",
                "test_audio.pcm",
                slice_u8.len()
            );
        }
        Err(e) => {
            println!(
                "Failed writing Audio file to storage at {}, reason: {}",
                "test_audio.cpm", e
            );
        }
    }
}
