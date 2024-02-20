use crate::cic::CicFilter;

use crate::fir::FloatFirDecimate;

/// Hard-coded FIR filter, designed for 96k sampling freq, with 8k cut-off, for
/// pre-decimation filter on last round of decimation
const LOWPASS_COEFFS: [f32; 50] = [
    2.25311824e-04,
    5.28281130e-03,
    5.56600120e-03,
    6.02832048e-03,
    4.17572036e-03,
    8.07054871e-05,
    -5.13611341e-03,
    -9.41126947e-03,
    -1.04681413e-02,
    -6.85823472e-03,
    1.06846141e-03,
    1.08391130e-02,
    1.83727875e-02,
    1.94600250e-02,
    1.16355776e-02,
    -4.21690294e-03,
    -2.32930647e-02,
    -3.78186694e-02,
    -3.93643088e-02,
    -2.18421132e-02,
    1.58325913e-02,
    6.85293828e-02,
    1.25715174e-01,
    1.74201321e-01,
    2.02001151e-01,
    2.02001151e-01,
    1.74201321e-01,
    1.25715174e-01,
    6.85293828e-02,
    1.58325913e-02,
    -2.18421132e-02,
    -3.93643088e-02,
    -3.78186694e-02,
    -2.32930647e-02,
    -4.21690294e-03,
    1.16355776e-02,
    1.94600250e-02,
    1.83727875e-02,
    1.08391130e-02,
    1.06846141e-03,
    -6.85823472e-03,
    -1.04681413e-02,
    -9.41126947e-03,
    -5.13611341e-03,
    8.07054871e-05,
    4.17572036e-03,
    6.02832048e-03,
    5.56600120e-03,
    5.28281130e-03,
    2.25311824e-04,
];

/// PdmProcessing struct stores filter state for decimation of all channels
const DEC1: usize = 8; // First stage decimation ratio
const DEC2: usize = 4; // Second stage decimation ratio
const DEC3: usize = 4; // Third stage decimation ratio
const ORDER1: usize = 4; // Order of first stage CIC
const ORDER2: usize = 3; // Order of second stage CIC

const SAMPLE_RATIO: usize = DEC1 * DEC2 * DEC3;

/// Number of samples at a time to pass to FIR filter Larger batch increases memory consumption, but
/// gives better performance up to a point
const BATCH_SIZE: usize = 64;

pub trait PdmProcessor {
    fn process_pdm(&self, pdm: &[u8], channel: usize, out: &mut [f32]);
}

pub struct PdmFilters {
    cic1: CicFilter<DEC1, ORDER1>,
    cic2: CicFilter<DEC2, ORDER2>,
    fir: FloatFirDecimate<{ LOWPASS_COEFFS.len() }, BATCH_SIZE, DEC3>,
    fir_buffer: [f32; BATCH_SIZE],
}

impl PdmFilters {
    pub fn new() -> Self {
        Self {
            cic1: CicFilter::new(),
            cic2: CicFilter::new(),
            fir: FloatFirDecimate::new(LOWPASS_COEFFS),
            fir_buffer: [0.0; BATCH_SIZE],
        }
    }

    pub fn process_pdm(&mut self, pdm: &[u8], channel: usize, out: &mut [f32]) {
        println!("Got pdm {} and out {}", pdm.len(), out.len());
        let total_samples = pdm.len() * 8 / 1 / SAMPLE_RATIO;
        let mut skip_samples = 0;
        let mut outpos: usize = 0;
        let mut batchpos: usize = 0;
        //  Deref to get around the RefMut and get a real reference
        // https://stackoverflow.com/questions/47060266/error-while-trying-to-borrow-2-fields-from-a-struct-wrapped-in-refcell
        let mut cic1 = self.cic1;
        let mut cic2 = self.cic2;
        let mut fir = self.fir;
        let mut fir_buffer = self.fir_buffer;
        cic1.process_pdm_buffer::<_, 1>(channel, pdm, |sample1| {
            cic2.push_sample(sample1, |sample2| {
                // Convert to float
                // Scale so that full-scale input results in +/- 1.0 output
                const FULL_SCALE: usize =
                    usize::pow(DEC1, ORDER1 as u32) * usize::pow(DEC2, ORDER2 as u32);
                let float_sample = sample2 as f32 / FULL_SCALE as f32;

                fir_buffer[batchpos] = float_sample;
                batchpos += 1;
                if batchpos == fir_buffer.len() {
                    // Low pass
                    fir.process_block(
                        fir_buffer.as_mut(),
                        &mut out[outpos..outpos + BATCH_SIZE / DEC3],
                    );
                    outpos += BATCH_SIZE / DEC3;

                    if skip_samples > 0 {
                        if skip_samples >= BATCH_SIZE / DEC3 {
                            // just skip the whole batch
                            outpos = 0;
                            skip_samples -= BATCH_SIZE / DEC3;
                        } else {
                            // We have to shift the data
                            for i in 0..BATCH_SIZE / DEC3 - skip_samples {
                                out[i] = out[i + skip_samples];
                            }
                            outpos -= skip_samples;
                            skip_samples = 0;
                        }
                    }
                    batchpos = 0;
                }
            });
        });

        // Finish any remaining partial batch
        if batchpos > 0 {
            fir.process_block(
                &fir_buffer[0..batchpos],
                &mut out[outpos..outpos + batchpos / DEC3],
            );
        }
        println!("Last pos {}", outpos);
    }
}
