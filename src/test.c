// Hello world! Cplayground is an online sandbox that makes it easy to try out
// code.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int DECIMATION_MAX = 128;
int SINCN = 3;


void convolve(uint32_t Signal[/* SignalLen */], unsigned short SignalLen,
              uint32_t Kernel[/* KernelLen */], unsigned short KernelLen,
              uint32_t Result[/* SignalLen + KernelLen - 1 */])
{
  uint16_t n;
 
  for (n = 0; n < SignalLen + KernelLen - 1; n++)
  {
    unsigned short kmin, kmax, k;
    
    Result[n] = 0;
    
    kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
    kmax = (n < SignalLen - 1) ? n : SignalLen - 1;
    
    for (k = kmin; k <= kmax; k++) {
      Result[n] += Signal[k] * Kernel[n - k];
    }
  }
}
int main() {
  uint32_t sinc[DECIMATION_MAX * SINCN];
    int decimation = 64;

  uint32_t sinc1[DECIMATION_MAX];
  uint32_t sinc2[DECIMATION_MAX * 2];
   for (int i = 0; i < decimation; i++) {
    sinc1[i] = 1;
  }
  sinc[0] = 0;
   int64_t sum = 0;
uint32_t coef[SINCN][DECIMATION_MAX];
int64_t sub_const = 0;

  sinc[decimation * SINCN - 1] = 0;      
  convolve(sinc1, decimation, sinc1, decimation, sinc2);
  convolve(sinc2, decimation * 2 - 1, sinc1, decimation, &sinc[1]);  
//   for(int i = 0; i < 200; i++){
//       printf("%d = %d\n",i+1,sinc[i+1]);
//   }
    int index = 0;
    for(int j = 0; j < SINCN; j++) {
    for (int i = 0; i < decimation; i++) {
      coef[j][i] = sinc[j * decimation + i];
      sum += sinc[j * decimation + i];

      printf("%d = %d\n",index,sinc[j * decimation + i]);             
      index +=1;

    }
    // break;
  }
 int32_t lut[256][DECIMATION_MAX / 8][SINCN];

  sub_const = sum >> 1;
  int FILTER_GAIN    = 16;
  int MaxVolume = 64;
  int div_const = sub_const * MaxVolume / 32768 / FILTER_GAIN;
  div_const = (div_const == 0 ? 1 : div_const);
  printf("Sum %d Sub const %d div %d\n",sum,sub_const,div_const);
  
    uint16_t c, d, s;
  for (s = 0; s < SINCN; s++)
  {
    uint32_t *coef_p = &coef[s][0];
    for (c = 0; c < 256; c++){
      for (d = 0; d < decimation / 8; d++){
        lut[c][d][s] = ((c >> 7)       ) * coef_p[d * 8    ] +
                       ((c >> 6) & 0x01) * coef_p[d * 8 + 1] +
                       ((c >> 5) & 0x01) * coef_p[d * 8 + 2] +
                       ((c >> 4) & 0x01) * coef_p[d * 8 + 3] +
                       ((c >> 3) & 0x01) * coef_p[d * 8 + 4] +
                       ((c >> 2) & 0x01) * coef_p[d * 8 + 5] +
                       ((c >> 1) & 0x01) * coef_p[d * 8 + 6] +
                       ((c     ) & 0x01) * coef_p[d * 8 + 7];
       printf("s = %d c=%d d= %d  %d\n",s,c,d,lut[c][d][s]);
          
      }
    }
  }
  
  return 0;
}


