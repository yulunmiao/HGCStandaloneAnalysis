#ifndef HGCSSDetector_h
#define HGCSSDetector_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TH2D.h"

enum DetectorEnum {
  FECAL,
  MECAL,
  BECAL,
  FHCAL,
  BHCAL1,
  BHCAL2
};

class HGCSSSubDetector {

public:
  HGCSSSubDetector():
    type(FECAL),
    name(""),
    layerIdMin(0),
    layerIdMax(0),
    mipWeight(1),
    absWeight(1),
    gevWeight(1),
    gevOffset(0),
    isSi(false),
    isScint(false),
    radiusLim(0)
  {};
  virtual ~HGCSSSubDetector(){};

  DetectorEnum type;
  std::string name;
  unsigned layerIdMin;
  unsigned layerIdMax;
  double mipWeight;
  double absWeight;
  double gevWeight;
  double gevOffset;
  bool isSi;
  bool isScint;
  double radiusLim;

  inline unsigned nLayers() const{
    return (layerIdMax-layerIdMin);
  };

private:

};

class HGCSSDetector {

public:
  friend HGCSSDetector & theDetector();

  inline void initialiseIndices(const unsigned versionNumber, const unsigned model=2){
    
    indices_.clear();
    indices_.resize(7,0);
    //fill layer indices
    if (versionNumber==22){
      indices_[4] = 0;
      indices_[5] = 10;
      indices_[6] = 10;
    }
    else if (versionNumber==28 || versionNumber==32) {
      indices_[4] = 0;
      indices_[5] = 12;
      indices_[6] = 12;
    }
    else if (versionNumber==23) {
      indices_[3] = 0;
      indices_[4] = 38;
      indices_[5] = 47;
      indices_[6] = 54;
    }
    else if (versionNumber==21) {
      indices_[3] = 0;
      indices_[4] = 24;
      indices_[5] = 34;
      indices_[6] = 34;
    }
    else if (versionNumber==27 || versionNumber==31) {
      indices_[3] = 0;
      indices_[4] = 12;
      indices_[5] = 24;
      indices_[6] = 24;
    }
    else if (versionNumber==38) {
      indices_[3] = 0;
      indices_[4] = 11;
      indices_[5] = 23;
      indices_[6] = 23;
    }
    else if (versionNumber==39) {
      indices_[3] = 0;
      indices_[4] = 9;
      indices_[5] = 21;
      indices_[6] = 21;
    }
    else if (versionNumber < 20){
      indices_[0] = 0;
      indices_[1] = versionNumber==8?11:10;
      indices_[2] = versionNumber==8?21:20;
      indices_[3] = versionNumber==8?31:30;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 30 || versionNumber == 60  || versionNumber==68 || (versionNumber >= 100 && versionNumber < 104)){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];

      sensitiveZ_.resize(indices_[3],0);
      sensitiveZ_[0] = 3198;
      sensitiveZ_[1] = 3207.1;
      sensitiveZ_[2] = 3222.4;
      sensitiveZ_[3] = 3231.5;
      sensitiveZ_[4] = 3246.8;
      sensitiveZ_[5] = 3255.9;
      sensitiveZ_[6] = 3271.2;
      sensitiveZ_[7] = 3280.3;
      sensitiveZ_[8] = 3295.6;
      sensitiveZ_[9] = 3304.7;
      sensitiveZ_[10] = 3320;
      sensitiveZ_[11] = 3329.1;
      sensitiveZ_[12] = 3344.4;
      sensitiveZ_[13] = 3353.5;
      sensitiveZ_[14] = 3368.8;
      sensitiveZ_[15] = 3377.9;
      sensitiveZ_[16] = 3393.2;
      sensitiveZ_[17] = 3402.3;
      sensitiveZ_[18] = 3417.6;
      sensitiveZ_[19] = 3426.7;
      sensitiveZ_[20] = 3442;
      sensitiveZ_[21] = 3451.1;
      sensitiveZ_[22] = 3466.4;
      sensitiveZ_[23] = 3475.5;
      sensitiveZ_[24] = 3490.8;
      sensitiveZ_[25] = 3499.9;
      sensitiveZ_[26] = 3515.2;
      sensitiveZ_[27] = 3524.3;
      if (model==3){
	sensitiveZ_[0] = -77.3;
	sensitiveZ_[1] = -68.2;
	sensitiveZ_[2] = -52.9;
	sensitiveZ_[3] = -43.8;
	sensitiveZ_[4] = -28.5;
	sensitiveZ_[5] = -19.4;
	sensitiveZ_[6] = -4.1;
	sensitiveZ_[7] = 5;
	sensitiveZ_[8] = 20.3;
	sensitiveZ_[9] = 29.4;
	sensitiveZ_[10] = 44.7;
	sensitiveZ_[11] = 53.8;
	sensitiveZ_[12] = 69.1;
	sensitiveZ_[13] = 78.2;
	sensitiveZ_[14] = 93.5;
	sensitiveZ_[15] = 102.6;
	sensitiveZ_[16] = 117.9;
	sensitiveZ_[17] = 127;
	sensitiveZ_[18] = 142.3;
	sensitiveZ_[19] = 151.4;
	sensitiveZ_[20] = 166.7;
	sensitiveZ_[21] = 175.8;
	sensitiveZ_[22] = 191.1;
	sensitiveZ_[23] = 200.2;
	sensitiveZ_[24] = 215.5;
	sensitiveZ_[25] = 224.6;
	sensitiveZ_[26] = 239.9;
	sensitiveZ_[27] = 249;
      }
    } 
    else if (versionNumber == 64){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];

      sensitiveZ_.resize(indices_[3],0);
      sensitiveZ_[0] = -96.05;
      sensitiveZ_[1] = -86.95;
      sensitiveZ_[2] = -68.65;
      sensitiveZ_[3] = -59.55;
      sensitiveZ_[4] = -41.25;
      sensitiveZ_[5] = -32.15;
      sensitiveZ_[6] = -13.85;
      sensitiveZ_[7] = -4.75;
      sensitiveZ_[8] = 13.55;
      sensitiveZ_[9] = 22.65;
      sensitiveZ_[10] = 40.95;
      sensitiveZ_[11] = 50.05;
      sensitiveZ_[12] = 68.35;
      sensitiveZ_[13] = 77.45;
      sensitiveZ_[14] = 95.75;
      sensitiveZ_[15] = 104.85;
      sensitiveZ_[16] = 123.15;
      sensitiveZ_[17] = 132.25;
      sensitiveZ_[18] = 150.55;
      sensitiveZ_[19] = 159.65;
      sensitiveZ_[20] = 177.95;
      sensitiveZ_[21] = 187.05;
      sensitiveZ_[22] = 205.35;
      sensitiveZ_[23] = 214.45;
      sensitiveZ_[24] = 232.75;
      sensitiveZ_[25] = 241.85;
      sensitiveZ_[26] = 260.15;
      sensitiveZ_[27] = 269.25;
    }
    else if (versionNumber == 67){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];

      sensitiveZ_.resize(indices_[3],0);
      sensitiveZ_[0] = 3200.5;
      sensitiveZ_[1] = 3209.6;
      sensitiveZ_[2] = 3229.9;
      sensitiveZ_[3] = 3239;
      sensitiveZ_[4] = 3259.3;
      sensitiveZ_[5] = 3268.4;
      sensitiveZ_[6] = 3288.7;
      sensitiveZ_[7] = 3297.8;
      sensitiveZ_[8] = 3318.1;
      sensitiveZ_[9] = 3327.2;
      sensitiveZ_[10] = 3347.5;
      sensitiveZ_[11] = 3356.6;
      sensitiveZ_[12] = 3376.9;
      sensitiveZ_[13] = 3386;
      sensitiveZ_[14] = 3406.3;
      sensitiveZ_[15] = 3415.4;
      sensitiveZ_[16] = 3435.7;
      sensitiveZ_[17] = 3444.8;
      sensitiveZ_[18] = 3465.1;
      sensitiveZ_[19] = 3474.2;
      sensitiveZ_[20] = 3494.5;
      sensitiveZ_[21] = 3503.6;
      sensitiveZ_[22] = 3523.9;
      sensitiveZ_[23] = 3533;
      sensitiveZ_[24] = 3553.3;
      sensitiveZ_[25] = 3562.4;
      sensitiveZ_[26] = 3582.7;
      sensitiveZ_[27] = 3591.8;
      if (model==3){
	sensitiveZ_[0] = -108.55;
	sensitiveZ_[1] = -99.45;
	sensitiveZ_[2] = -79.15;
	sensitiveZ_[3] = -70.05;
	sensitiveZ_[4] = -49.75;
	sensitiveZ_[5] = -40.65;
	sensitiveZ_[6] = -20.35;
	sensitiveZ_[7] = -11.25;
	sensitiveZ_[8] = 9.05;
	sensitiveZ_[9] = 18.15;
	sensitiveZ_[10] = 38.45;
	sensitiveZ_[11] = 47.55;
	sensitiveZ_[12] = 67.85;
	sensitiveZ_[13] = 76.95;
	sensitiveZ_[14] = 97.25;
	sensitiveZ_[15] = 106.35;
	sensitiveZ_[16] = 126.65;
	sensitiveZ_[17] = 135.75;
	sensitiveZ_[18] = 156.05;
	sensitiveZ_[19] = 165.15;
	sensitiveZ_[20] = 185.45;
	sensitiveZ_[21] = 194.55;
	sensitiveZ_[22] = 214.85;
	sensitiveZ_[23] = 223.95;
	sensitiveZ_[24] = 244.25;
	sensitiveZ_[25] = 253.35;
	sensitiveZ_[26] = 273.65;
	sensitiveZ_[27] = 282.75;
      }
    }
    else if (versionNumber == 65){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];

      sensitiveZ_.resize(indices_[3],0);
      sensitiveZ_[0] = -78.7;
      sensitiveZ_[1] = -70.4;
      sensitiveZ_[2] = -50.3;
      sensitiveZ_[3] = -42;
      sensitiveZ_[4] = -21.9;
      sensitiveZ_[5] = -13.6;
      sensitiveZ_[6] = 6.5;
      sensitiveZ_[7] = 14.8;
      sensitiveZ_[8] = 34.9;
      sensitiveZ_[9] = 43.2;
      sensitiveZ_[10] = 63.3;
      sensitiveZ_[11] = 71.6;
      sensitiveZ_[12] = 91.7;
      sensitiveZ_[13] = 100;
      sensitiveZ_[14] = 120.1;
      sensitiveZ_[15] = 128.4;
      sensitiveZ_[16] = 148.5;
      sensitiveZ_[17] = 156.8;
      sensitiveZ_[18] = 176.9;
      sensitiveZ_[19] = 185.2;
      sensitiveZ_[20] = 205.3;
      sensitiveZ_[21] = 213.6;
      sensitiveZ_[22] = 233.7;
      sensitiveZ_[23] = 242;
      sensitiveZ_[24] = 262.1;
      sensitiveZ_[25] = 270.4;
      sensitiveZ_[26] = 290.5;
      sensitiveZ_[27] = 298.8;
     }
     else if (versionNumber == 66){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 24;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];

      sensitiveZ_.resize(indices_[3],0);
      sensitiveZ_[0] = -59.3;
      sensitiveZ_[1] = -50.2;
      sensitiveZ_[2] = -29.3;
      sensitiveZ_[3] = -20.2;
      sensitiveZ_[4] = 0.7;
      sensitiveZ_[5] = 9.8;
      sensitiveZ_[6] = 30.7;
      sensitiveZ_[7] = 39.8;
      sensitiveZ_[8] = 60.7;
      sensitiveZ_[9] = 69.8;
      sensitiveZ_[10] = 90.7;
      sensitiveZ_[11] = 99.8;
      sensitiveZ_[12] = 120.7;
      sensitiveZ_[13] = 129.8;
      sensitiveZ_[14] = 150.7;
      sensitiveZ_[15] = 159.8;
      sensitiveZ_[16] = 180.7;
      sensitiveZ_[17] = 189.8;
      sensitiveZ_[18] = 210.7;
      sensitiveZ_[19] = 219.8;
      sensitiveZ_[20] = 240.7;
      sensitiveZ_[21] = 249.8;
      sensitiveZ_[22] = 270.7;
      sensitiveZ_[23] = 279.8;
     }
     else if (versionNumber == 70) {
       indices_[0] = 0;
       indices_[1] = 10;
       indices_[2] = 18;
       indices_[3] = 26;
       indices_[4] = indices_[3];
       indices_[5] = indices_[3];
       indices_[6] = indices_[3];
       
       sensitiveZ_.resize(indices_[3],0);
       sensitiveZ_[0] = 3213.95;
       sensitiveZ_[1] = 3223.7;
       sensitiveZ_[2] = 3244.9;
       sensitiveZ_[3] = 3254.65;
       sensitiveZ_[4] = 3275.85;
       sensitiveZ_[5] = 3285.6;
       sensitiveZ_[6] = 3306.8;
       sensitiveZ_[7] = 3316.55;
       sensitiveZ_[8] = 3337.75;
       sensitiveZ_[9] = 3347.5;
       sensitiveZ_[10] = 3368.7;
       sensitiveZ_[11] = 3378.45;
       sensitiveZ_[12] = 3399.65;
       sensitiveZ_[13] = 3409.4;
       sensitiveZ_[14] = 3430.6;
       sensitiveZ_[15] = 3440.35;
       sensitiveZ_[16] = 3461.55;
       sensitiveZ_[17] = 3471.3;
       sensitiveZ_[18] = 3495.52;
       sensitiveZ_[19] = 3505.27;
       sensitiveZ_[20] = 3529.49;
       sensitiveZ_[21] = 3539.24;
       sensitiveZ_[22] = 3563.46;
       sensitiveZ_[23] = 3573.21;
       sensitiveZ_[24] = 3597.43;
       sensitiveZ_[25] = 3607.18;

       if (model==3){
	 sensitiveZ_[0] = -106.465;
	 sensitiveZ_[1] = -96.715;
	 sensitiveZ_[2] = -75.515;
	 sensitiveZ_[3] = -65.765;
	 sensitiveZ_[4] = -44.565;
	 sensitiveZ_[5] = -34.815;
	 sensitiveZ_[6] = -13.615;
	 sensitiveZ_[7] = -3.865;
	 sensitiveZ_[8] = 17.335;
	 sensitiveZ_[9] = 27.085;
	 sensitiveZ_[10] = 48.285;
	 sensitiveZ_[11] = 58.035;
	 sensitiveZ_[12] = 79.235;
	 sensitiveZ_[13] = 88.985;
	 sensitiveZ_[14] = 110.185;
	 sensitiveZ_[15] = 119.935;
	 sensitiveZ_[16] = 141.135;
	 sensitiveZ_[17] = 150.885;
	 sensitiveZ_[18] = 175.105;
	 sensitiveZ_[19] = 184.855;
	 sensitiveZ_[20] = 209.075;
	 sensitiveZ_[21] = 218.825;
	 sensitiveZ_[22] = 243.045;
	 sensitiveZ_[23] = 252.795;
	 sensitiveZ_[24] = 277.015;
	 sensitiveZ_[25] = 286.765;
       }
       
       etaBoundary_.resize(indices_[3],0);
       etaBoundary_[0] =1.49186;
       etaBoundary_[1] =1.49327;
       etaBoundary_[2] =1.49467;
       etaBoundary_[3] =1.49605;
       etaBoundary_[4] =1.49743;
       etaBoundary_[5] =1.49879;
       etaBoundary_[6] =1.50015;
       etaBoundary_[7] =1.50149;
       etaBoundary_[8] =1.50283;
       etaBoundary_[9] =1.50415;
       etaBoundary_[10] =1.50547;
       etaBoundary_[11] =1.50677;
       etaBoundary_[12] =1.50808;
       etaBoundary_[13] =1.50935;
       etaBoundary_[14] =1.51064;
       etaBoundary_[15] =1.5119;
       etaBoundary_[16] =1.51316;
       etaBoundary_[17] =1.51441;
       etaBoundary_[18] =1.51589;
       etaBoundary_[19] =1.51711;
       etaBoundary_[20] =1.51858;
       etaBoundary_[21] =1.51978;
       etaBoundary_[22] =1.52122;
       etaBoundary_[23] =1.52241;
       etaBoundary_[24] =1.52382;
       etaBoundary_[25] =1.52843;
     }
     else if (versionNumber == 33){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = 40;
      indices_[5] = 52;
      indices_[6] = 52;
    }
    else if (versionNumber == 34){
      indices_[0] = 0;
      indices_[1] = 8;
      indices_[2] = 16;
      indices_[3] = 24;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 35){
      indices_[0] = 0;
      indices_[1] = 6;
      indices_[2] = 12;
      indices_[3] = 18;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 36){
      indices_[0] = 0;
      indices_[1] = 8;
      indices_[2] = 16;
      indices_[3] = 24;
      indices_[4] = 35;
      indices_[5] = 47;
      indices_[6] = 47;
    }
    else if (versionNumber == 37){
      indices_[0] = 0;
      indices_[1] = 6;
      indices_[2] = 12;
      indices_[3] = 18;
      indices_[4] = 27;
      indices_[5] = 39;
      indices_[6] = 39;
    }
    else if (versionNumber == 61){
      indices_[3] = 0;
      indices_[4] = 25;
      indices_[5] = 41;
      indices_[6] = 41;

      sensitiveZ_.resize(indices_[6],0);
      sensitiveZ_[0] = 3037;
      sensitiveZ_[1] = 3086;
      sensitiveZ_[2] = 3135;
      sensitiveZ_[3] = 3184;
      sensitiveZ_[4] = 3233;
      sensitiveZ_[5] = 3282;
      sensitiveZ_[6] = 3331;
      sensitiveZ_[7] = 3380;
      sensitiveZ_[8] = 3429;
      sensitiveZ_[9] = 3479.9;
      sensitiveZ_[10] = 3530.8;
      sensitiveZ_[11] = 3581.7;
      sensitiveZ_[12] = 3665.6;
      sensitiveZ_[13] = 3749.5;
      sensitiveZ_[14] = 3833.4;
      sensitiveZ_[15] = 3917.3;
      sensitiveZ_[16] = 4001.2;
      sensitiveZ_[17] = 4085.1;
      sensitiveZ_[18] = 4169;
      sensitiveZ_[19] = 4252.9;
      sensitiveZ_[20] = 4336.8;
      sensitiveZ_[21] = 4420.7;
      sensitiveZ_[22] = 4504.6;
      sensitiveZ_[23] = 4588.5;
      sensitiveZ_[24] = 0;
      sensitiveZ_[25] = 3430.8;
      sensitiveZ_[26] = 3481.7;
      sensitiveZ_[27] = 3532.6;
      sensitiveZ_[28] = 3583.5;
      sensitiveZ_[29] = 3667.4;
      sensitiveZ_[30] = 3751.3;
      sensitiveZ_[31] = 3835.2;
      sensitiveZ_[32] = 3919.1;
      sensitiveZ_[33] = 4003;
      sensitiveZ_[34] = 4086.9;
      sensitiveZ_[35] = 4170.8;
      sensitiveZ_[36] = 4254.7;
      sensitiveZ_[37] = 4338.6;
      sensitiveZ_[38] = 4422.5;
      sensitiveZ_[39] = 4506.4;
      sensitiveZ_[40] = 4590.3;
      etaBoundary_.resize(indices_[6],0);
      for (unsigned iL(0); iL<8; ++iL){
	etaBoundary_[iL] = 1.4;
      }
      etaBoundary_[8] = 1.72042;
      etaBoundary_[9] = 1.81718;
      etaBoundary_[10] = 1.82927;
      etaBoundary_[11] = 1.91612;
      etaBoundary_[12] = 2.02287;
      etaBoundary_[13] = 2.09617;
      etaBoundary_[14] = 2.21281;
      etaBoundary_[15] = 2.28463;
      etaBoundary_[16] = 2.30324;
      etaBoundary_[17] = 2.32153;
      etaBoundary_[18] = 2.33949;
      etaBoundary_[19] = 2.35714;
      etaBoundary_[20] = 2.37449;
      etaBoundary_[21] = 2.39155;
      etaBoundary_[22] = 2.40832;
      etaBoundary_[23] = 2.42483;

      etaBoundary_[25] = 1.72042;
      etaBoundary_[26] = 1.81718;
      etaBoundary_[27] = 1.82927;
      etaBoundary_[28] = 1.91612;
      etaBoundary_[29] = 2.02287;
      etaBoundary_[30] = 2.09617;
      etaBoundary_[31] = 2.21281;
      etaBoundary_[32] = 2.28463;
      etaBoundary_[33] = 2.30324;
      etaBoundary_[34] = 2.32153;
      etaBoundary_[35] = 2.33949;
      etaBoundary_[36] = 2.35714;
      etaBoundary_[37] = 2.37449;
      etaBoundary_[38] = 2.39155;
      etaBoundary_[39] = 2.40832;
      etaBoundary_[40] = 2.42483;
    }
    else if (versionNumber == 62){
      indices_[4] = 0;
      indices_[5] = 16;
      indices_[6] = 16;

      sensitiveZ_.resize(indices_[6],0);
      sensitiveZ_[0] = 3040.5;
      sensitiveZ_[1] = 3091.4;
      sensitiveZ_[2] = 3142.3;
      sensitiveZ_[3] = 3193.2;
      sensitiveZ_[4] = 3277.1;
      sensitiveZ_[5] = 3361;
      sensitiveZ_[6] = 3444.9;
      sensitiveZ_[7] = 3528.8;
      sensitiveZ_[8] = 3612.7;
      sensitiveZ_[9] = 3696.6;
      sensitiveZ_[10] = 3780.5;
      sensitiveZ_[11] = 3864.4;
      sensitiveZ_[12] = 3948.3;
      sensitiveZ_[13] = 4032.2;
      sensitiveZ_[14] = 4116.1;
      sensitiveZ_[15] = 4200;

      etaBoundary_.resize(indices_[6],0);
      etaBoundary_[0] = 1.72042;
      etaBoundary_[1] = 1.81718;
      etaBoundary_[2] = 1.82927;
      etaBoundary_[3] = 1.91612;
      etaBoundary_[4] = 2.02287;
      etaBoundary_[5] = 2.09617;
      etaBoundary_[6] = 2.21281;
      etaBoundary_[7] = 2.28463;
      etaBoundary_[8] = 2.30324;
      etaBoundary_[9] = 2.32153;
      etaBoundary_[10] = 2.33949;
      etaBoundary_[11] = 2.35714;
      etaBoundary_[12] = 2.37449;
      etaBoundary_[13] = 2.39155;
      etaBoundary_[14] = 2.40832;
      etaBoundary_[15] = 2.42483;
    }
    else if (versionNumber == 63){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = 53;
      indices_[5] = 57;
      indices_[6] = 69;
      sensitiveZ_.resize(indices_[6],0);
      sensitiveZ_[0] = 3198;
      sensitiveZ_[1] = 3207.1;
      sensitiveZ_[2] = 3222.4;
      sensitiveZ_[3] = 3231.5;
      sensitiveZ_[4] = 3246.8;
      sensitiveZ_[5] = 3255.9;
      sensitiveZ_[6] = 3271.2;
      sensitiveZ_[7] = 3280.3;
      sensitiveZ_[8] = 3295.6;
      sensitiveZ_[9] = 3304.7;
      sensitiveZ_[10] = 3320;
      sensitiveZ_[11] = 3329.1;
      sensitiveZ_[12] = 3344.4;
      sensitiveZ_[13] = 3353.5;
      sensitiveZ_[14] = 3368.8;
      sensitiveZ_[15] = 3377.9;
      sensitiveZ_[16] = 3393.2;
      sensitiveZ_[17] = 3402.3;
      sensitiveZ_[18] = 3417.6;
      sensitiveZ_[19] = 3426.7;
      sensitiveZ_[20] = 3442;
      sensitiveZ_[21] = 3451.1;
      sensitiveZ_[22] = 3466.4;
      sensitiveZ_[23] = 3475.5;
      sensitiveZ_[24] = 3490.8;
      sensitiveZ_[25] = 3499.9;
      sensitiveZ_[26] = 3515.2;
      sensitiveZ_[27] = 3524.3;
      sensitiveZ_[28] = 3577.4;
      sensitiveZ_[29] = 3626.4;
      sensitiveZ_[30] = 3675.4;
      sensitiveZ_[31] = 3724.4;
      sensitiveZ_[32] = 3773.4;
      sensitiveZ_[33] = 3822.4;
      sensitiveZ_[34] = 3871.4;
      sensitiveZ_[35] = 3920.4;
      sensitiveZ_[36] = 3969.4;
      sensitiveZ_[37] = 4020.3;
      sensitiveZ_[38] = 4071.2;
      sensitiveZ_[39] = 4122.1;
      sensitiveZ_[40] = 4206;
      sensitiveZ_[41] = 4289.9;
      sensitiveZ_[42] = 4373.8;
      sensitiveZ_[43] = 4457.7;
      sensitiveZ_[44] = 4541.6;
      sensitiveZ_[45] = 4625.5;
      sensitiveZ_[46] = 4709.4;
      sensitiveZ_[47] = 4793.3;
      sensitiveZ_[48] = 4877.2;
      sensitiveZ_[49] = 4961.1;
      sensitiveZ_[50] = 5045;
      sensitiveZ_[51] = 5128.9;
      sensitiveZ_[52] = 0;
      sensitiveZ_[53] = 3971.2;
      sensitiveZ_[54] = 4022.1;
      sensitiveZ_[55] = 4073;
      sensitiveZ_[56] = 4123.9;
      sensitiveZ_[57] = 4207.8;
      sensitiveZ_[58] = 4291.7;
      sensitiveZ_[59] = 4375.6;
      sensitiveZ_[60] = 4459.5;
      sensitiveZ_[61] = 4543.4;
      sensitiveZ_[62] = 4627.3;
      sensitiveZ_[63] = 4711.2;
      sensitiveZ_[64] = 4795.1;
      sensitiveZ_[65] = 4879;
      sensitiveZ_[66] = 4962.9;
      sensitiveZ_[67] = 5046.8;
      sensitiveZ_[68] = 5130.7;

      etaBoundary_.resize(indices_[6],0);
      etaBoundary_[0]  = 1.461;
      etaBoundary_[1]  = 1.464;
      etaBoundary_[2]  = 1.463;
      etaBoundary_[3]  = 1.466;
      etaBoundary_[4]  = 1.465;
      etaBoundary_[5]  = 1.468;
      etaBoundary_[6]  = 1.467;
      etaBoundary_[7]  = 1.469;
      etaBoundary_[8]  = 1.469;
      etaBoundary_[9]  = 1.471;
      etaBoundary_[10]  = 1.471;
      etaBoundary_[11]  = 1.473;
      etaBoundary_[12]  = 1.472;
      etaBoundary_[13]  = 1.475;
      etaBoundary_[14]  = 1.474;
      etaBoundary_[15]  = 1.477;
      etaBoundary_[16]  = 1.476;
      etaBoundary_[17]  = 1.478;
      etaBoundary_[18]  = 1.478;
      etaBoundary_[19]  = 1.480;
      etaBoundary_[20]  = 1.479;
      etaBoundary_[21]  = 1.482;
      etaBoundary_[22]  = 1.481;
      etaBoundary_[23]  = 1.483;
      etaBoundary_[24]  = 1.483;
      etaBoundary_[25]  = 1.485;
      etaBoundary_[26]  = 1.484;
      etaBoundary_[27]  = 1.487;
      etaBoundary_[28]  = 1.487;
      etaBoundary_[29]  = 1.490;
      etaBoundary_[30]  = 1.494;
      etaBoundary_[31]  = 1.497;
      etaBoundary_[32]  = 1.500;
      etaBoundary_[33]  = 1.503;
      etaBoundary_[34]  = 1.506;
      etaBoundary_[35]  = 1.493;
	
      etaBoundary_[36] = 1.72042;
      etaBoundary_[37] = 1.81718;
      etaBoundary_[38] = 1.82927;
      etaBoundary_[39] = 1.91612;
      etaBoundary_[40] = 2.02287;
      etaBoundary_[41] = 2.09617;
      etaBoundary_[42] = 2.21281;
      etaBoundary_[43] = 2.28463;
      etaBoundary_[44] = 2.30324;
      etaBoundary_[45] = 2.32153;
      etaBoundary_[46] = 2.33949;
      etaBoundary_[47] = 2.35714;
      etaBoundary_[48] = 2.37449;
      etaBoundary_[49] = 2.39155;
      etaBoundary_[50] = 2.40832;
      etaBoundary_[51] = 2.42483;

      etaBoundary_[53] = 1.72042;
      etaBoundary_[54] = 1.81718;
      etaBoundary_[55] = 1.82927;
      etaBoundary_[56] = 1.91612;
      etaBoundary_[57] = 2.02287;
      etaBoundary_[58] = 2.09617;
      etaBoundary_[59] = 2.21281;
      etaBoundary_[60] = 2.28463;
      etaBoundary_[61] = 2.30324;
      etaBoundary_[62] = 2.32153;
      etaBoundary_[63] = 2.33949;
      etaBoundary_[64] = 2.35714;
      etaBoundary_[65] = 2.37449;
      etaBoundary_[66] = 2.39155;
      etaBoundary_[67] = 2.40832;
      etaBoundary_[68] = 2.42483;

    }else if(versionNumber==69) {

      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = 53;
      indices_[5] = 57;
      indices_[6] = 69;

      sensitiveZ_.resize(indices_[6],0);
      sensitiveZ_[0] = 3198;
      sensitiveZ_[1] = 3207.1;
      sensitiveZ_[2] = 3222.4;
      sensitiveZ_[3] = 3231.5;
      sensitiveZ_[4] = 3246.8;
      sensitiveZ_[5] = 3255.9;
      sensitiveZ_[6] = 3271.2;
      sensitiveZ_[7] = 3280.3;
      sensitiveZ_[8] = 3295.6;
      sensitiveZ_[9] = 3304.7;
      sensitiveZ_[10] = 3320;
      sensitiveZ_[11] = 3329.1;
      sensitiveZ_[12] = 3344.4;
      sensitiveZ_[13] = 3353.5;
      sensitiveZ_[14] = 3368.8;
      sensitiveZ_[15] = 3377.9;
      sensitiveZ_[16] = 3393.2;
      sensitiveZ_[17] = 3402.3;
      sensitiveZ_[18] = 3417.6;
      sensitiveZ_[19] = 3426.7;
      sensitiveZ_[20] = 3442;
      sensitiveZ_[21] = 3451.1;
      sensitiveZ_[22] = 3466.4;
      sensitiveZ_[23] = 3475.5;
      sensitiveZ_[24] = 3490.8;
      sensitiveZ_[25] = 3499.9;
      sensitiveZ_[26] = 3515.2;
      sensitiveZ_[27] = 3524.3;
      sensitiveZ_[28] = 3577.4;
      sensitiveZ_[29] = 3626.4;
      sensitiveZ_[30] = 3675.4;
      sensitiveZ_[31] = 3724.4;
      sensitiveZ_[32] = 3773.4;
      sensitiveZ_[33] = 3822.4;
      sensitiveZ_[34] = 3871.4;
      sensitiveZ_[35] = 3920.4;
      sensitiveZ_[36] = 3969.4;
      sensitiveZ_[37] = 4020.3;
      sensitiveZ_[38] = 4071.2;
      sensitiveZ_[39] = 4122.1;
      sensitiveZ_[40] = 4206;
      sensitiveZ_[41] = 4289.9;
      sensitiveZ_[42] = 4373.8;
      sensitiveZ_[43] = 4457.7;
      sensitiveZ_[44] = 4541.6;
      sensitiveZ_[45] = 4625.5;
      sensitiveZ_[46] = 4709.4;
      sensitiveZ_[47] = 4793.3;
      sensitiveZ_[48] = 4877.2;
      sensitiveZ_[49] = 4961.1;
      sensitiveZ_[50] = 5045;
      sensitiveZ_[51] = 5128.9;
      sensitiveZ_[52] = 0;
      sensitiveZ_[53] = 3971.2;
      sensitiveZ_[54] = 4022.1;
      sensitiveZ_[55] = 4073;
      sensitiveZ_[56] = 4123.9;
      sensitiveZ_[57] = 4207.8;
      sensitiveZ_[58] = 4291.7;
      sensitiveZ_[59] = 4375.6;
      sensitiveZ_[60] = 4459.5;
      sensitiveZ_[61] = 4543.4;
      sensitiveZ_[62] = 4627.3;
      sensitiveZ_[63] = 4711.2;
      sensitiveZ_[64] = 4795.1;
      sensitiveZ_[65] = 4879;
      sensitiveZ_[66] = 4962.9;
      sensitiveZ_[67] = 5046.8;
      sensitiveZ_[68] = 5130.7;

      etaBoundary_.resize(indices_[6],0);
      etaBoundary_[0] =1.461;
      etaBoundary_[1] =1.464;
      etaBoundary_[2] =1.463;
      etaBoundary_[3] =1.466;
      etaBoundary_[4] =1.465;
      etaBoundary_[5] =1.468;
      etaBoundary_[6] =1.467;
      etaBoundary_[7] =1.469;
      etaBoundary_[8] =1.469;
      etaBoundary_[9] =1.471;
      etaBoundary_[10] =1.471;
      etaBoundary_[11] =1.473;
      etaBoundary_[12] =1.472;
      etaBoundary_[13] =1.475;
      etaBoundary_[14] =1.474;
      etaBoundary_[15] =1.477;
      etaBoundary_[16] =1.476;
      etaBoundary_[17] =1.478;
      etaBoundary_[18] =1.478;
      etaBoundary_[19] =1.48;
      etaBoundary_[20] =1.479;
      etaBoundary_[21] =1.482;
      etaBoundary_[22] =1.481;
      etaBoundary_[23] =1.483;
      etaBoundary_[24] =1.483;
      etaBoundary_[25] =1.485;
      etaBoundary_[26] =1.484;
      etaBoundary_[27] =1.487;
      etaBoundary_[28] =1.487;
      etaBoundary_[29] =1.49;
      etaBoundary_[30] =1.494;
      etaBoundary_[31] =1.497;
      etaBoundary_[32] =1.5;
      etaBoundary_[33] =1.503;
      etaBoundary_[34] =1.506;
      etaBoundary_[35] =1.493;
      etaBoundary_[36] =1.66415;
      etaBoundary_[37] =1.67572;
      etaBoundary_[38] =1.68761;
      etaBoundary_[39] =1.69936;
      etaBoundary_[40] =1.71099;
      etaBoundary_[41] =1.72989;
      etaBoundary_[42] =1.74846;
      etaBoundary_[43] =1.76671;
      etaBoundary_[44] =1.88338;
      etaBoundary_[45] =1.9012;
      etaBoundary_[46] =2.0303;
      etaBoundary_[47] =2.04768;
      etaBoundary_[48] =2.06477;
      etaBoundary_[49] =2.08158;
      etaBoundary_[50] =2.2441;
      etaBoundary_[51] =2.30169;
      etaBoundary_[52] =1.49193;
      etaBoundary_[53] =1.66415;
      etaBoundary_[54] =1.67572;
      etaBoundary_[55] =1.68761;
      etaBoundary_[56] =1.69936;
      etaBoundary_[57] =1.71099;
      etaBoundary_[58] =1.72989;
      etaBoundary_[59] =1.74846;
      etaBoundary_[60] =1.76671;
      etaBoundary_[61] =1.88338;
      etaBoundary_[62] =1.9012;
      etaBoundary_[63] =2.0303;
      etaBoundary_[64] =2.04768;
      etaBoundary_[65] =2.06477;
      etaBoundary_[66] =2.08158;
      etaBoundary_[67] =2.2441;
      etaBoundary_[68] =2.26051;
    }
    else if (versionNumber == 73){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 18;
      indices_[3] = 26;
      indices_[4] = 48;
      indices_[5] = 52;
      indices_[6] = 62;

      sensitiveZ_.resize(indices_[6],0);
      sensitiveZ_[0] = 3213.95;
      sensitiveZ_[1] = 3223.7;
      sensitiveZ_[2] = 3244.9;
      sensitiveZ_[3] = 3254.65;
      sensitiveZ_[4] = 3275.85;
      sensitiveZ_[5] = 3285.6;
      sensitiveZ_[6] = 3306.8;
      sensitiveZ_[7] = 3316.55;
      sensitiveZ_[8] = 3337.75;
      sensitiveZ_[9] = 3347.5;
      sensitiveZ_[10] = 3368.7;
      sensitiveZ_[11] = 3378.45;
      sensitiveZ_[12] = 3399.65;
      sensitiveZ_[13] = 3409.4;
      sensitiveZ_[14] = 3430.6;
      sensitiveZ_[15] = 3440.35;
      sensitiveZ_[16] = 3461.55;
      sensitiveZ_[17] = 3471.3;
      sensitiveZ_[18] = 3495.52;
      sensitiveZ_[19] = 3505.27;
      sensitiveZ_[20] = 3529.49;
      sensitiveZ_[21] = 3539.24;
      sensitiveZ_[22] = 3563.46;
      sensitiveZ_[23] = 3573.21;
      sensitiveZ_[24] = 3597.43;
      sensitiveZ_[25] = 3607.18;
      sensitiveZ_[26] = 3674.73;
      sensitiveZ_[27] = 3737.78;
      sensitiveZ_[28] = 3800.83;
      sensitiveZ_[29] = 3863.88;
      sensitiveZ_[30] = 3926.93;
      sensitiveZ_[31] = 3989.98;
      sensitiveZ_[32] = 4053.03;
      sensitiveZ_[33] = 4116.08;
      sensitiveZ_[34] = 4179.13;
      sensitiveZ_[35] = 4242.18;
      sensitiveZ_[36] = 4305.23;
      sensitiveZ_[37] = 4387.48;
      sensitiveZ_[38] = 4469.73;
      sensitiveZ_[39] = 4551.98;
      sensitiveZ_[40] = 4634.23;
      sensitiveZ_[41] = 4716.48;
      sensitiveZ_[42] = 4798.73;
      sensitiveZ_[43] = 4880.98;
      sensitiveZ_[44] = 4963.23;
      sensitiveZ_[45] = 5045.48;
      sensitiveZ_[46] = 5127.73;
      sensitiveZ_[47] = 0.;
      sensitiveZ_[48] = 4112.78;
      sensitiveZ_[49] = 4175.83;
      sensitiveZ_[50] = 4238.88;
      sensitiveZ_[51] = 4301.93;
      sensitiveZ_[52] = 4384.18;
      sensitiveZ_[53] = 4466.43;
      sensitiveZ_[54] = 4548.68;
      sensitiveZ_[55] = 4630.93;
      sensitiveZ_[56] = 4713.18;
      sensitiveZ_[57] = 4795.43;
      sensitiveZ_[58] = 4877.68;
      sensitiveZ_[59] = 4959.93;
      sensitiveZ_[60] = 5042.18;
      sensitiveZ_[61] = 5124.43;

      etaBoundary_.resize(indices_[6],0);
      etaBoundary_[0] =1.49186;
      etaBoundary_[1] =1.49327;
      etaBoundary_[2] =1.49467;
      etaBoundary_[3] =1.49605;
      etaBoundary_[4] =1.49743;
      etaBoundary_[5] =1.49879;
      etaBoundary_[6] =1.50015;
      etaBoundary_[7] =1.50149;
      etaBoundary_[8] =1.50283;
      etaBoundary_[9] =1.50415;
      etaBoundary_[10] =1.50547;
      etaBoundary_[11] =1.50677;
      etaBoundary_[12] =1.50808;
      etaBoundary_[13] =1.50935;
      etaBoundary_[14] =1.51064;
      etaBoundary_[15] =1.5119;
      etaBoundary_[16] =1.51316;
      etaBoundary_[17] =1.51441;
      etaBoundary_[18] =1.51589;
      etaBoundary_[19] =1.51711;
      etaBoundary_[20] =1.51858;
      etaBoundary_[21] =1.51978;
      etaBoundary_[22] =1.52122;
      etaBoundary_[23] =1.52241;
      etaBoundary_[24] =1.52382;
      etaBoundary_[25] =1.52499;
      etaBoundary_[26] =1.53001;
      etaBoundary_[27] =1.53457;
      etaBoundary_[28] =1.53899;
      etaBoundary_[29] =1.54269;
      etaBoundary_[30] =1.51664;
      etaBoundary_[31] =1.49222;
      etaBoundary_[32] =1.46929;
      etaBoundary_[33] =1.71129;
      etaBoundary_[34] =1.72552;
      etaBoundary_[35] =1.73956;
      etaBoundary_[36] =1.75342;
      etaBoundary_[37] =1.77123;
      etaBoundary_[38] =1.88752;
      etaBoundary_[39] =2.01639;
      etaBoundary_[40] =2.03366;
      etaBoundary_[41] =2.05064;
      etaBoundary_[42] =2.06735;
      etaBoundary_[43] =2.22966;
      etaBoundary_[44] =2.24597;
      etaBoundary_[45] =2.26202;
      etaBoundary_[46] =2.27783;
      etaBoundary_[47] =2.27783;
      etaBoundary_[48] =1.71129;
      etaBoundary_[49] =1.72552;
      etaBoundary_[50] =1.73956;
      etaBoundary_[51] =1.75342;
      etaBoundary_[52] =1.77123;
      etaBoundary_[53] =1.88752;
      etaBoundary_[54] =2.01639;
      etaBoundary_[55] =2.03366;
      etaBoundary_[56] =2.05064;
      etaBoundary_[57] =2.06735;
      etaBoundary_[58] =2.22966;
      etaBoundary_[59] =2.24597;
      etaBoundary_[60] =2.26202;
      etaBoundary_[61] =2.27783;

    }else if (versionNumber==80) {
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 26;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];

      sensitiveZ_.resize(indices_[3],0);
      sensitiveZ_[0] = 3213.95;
      sensitiveZ_[1] = 3223.94;
      sensitiveZ_[2] = 3245.64;
      sensitiveZ_[3] = 3255.64;
      sensitiveZ_[4] = 3277.34;
      sensitiveZ_[5] = 3287.33;
      sensitiveZ_[6] = 3309.03;
      sensitiveZ_[7] = 3319.03;
      sensitiveZ_[8] = 3340.73;
      sensitiveZ_[9] = 3350.72;
      sensitiveZ_[10] = 3372.42;
      sensitiveZ_[11] = 3382.42;
      sensitiveZ_[12] = 3404.12;
      sensitiveZ_[13] = 3414.11;
      sensitiveZ_[14] = 3435.81;
      sensitiveZ_[15] = 3445.8;
      sensitiveZ_[16] = 3467.5;
      sensitiveZ_[17] = 3477.5;
      sensitiveZ_[18] = 3499.2;
      sensitiveZ_[19] = 3509.19;
      sensitiveZ_[20] = 3530.89;
      sensitiveZ_[21] = 3540.89;
      sensitiveZ_[22] = 3562.59;
      sensitiveZ_[23] = 3572.58;
      sensitiveZ_[24] = 3594.28;
      sensitiveZ_[25] = 3604.28;
      if (model==3){
	sensitiveZ_[0] = -105.013;
	sensitiveZ_[1] = -95.019;
	sensitiveZ_[2] = -73.3185;
	sensitiveZ_[3] = -63.3246;
	sensitiveZ_[4] = -41.6242;
	sensitiveZ_[5] = -31.6303;
	sensitiveZ_[6] = -9.92988;
	sensitiveZ_[7] = 0.063983;
	sensitiveZ_[8] = 21.7644;
	sensitiveZ_[9] = 31.7583;
	sensitiveZ_[10] = 53.4588;
	sensitiveZ_[11] = 63.4526;
	sensitiveZ_[12] = 85.1531;
	sensitiveZ_[13] = 95.1469;
	sensitiveZ_[14] = 116.847;
	sensitiveZ_[15] = 126.841;
	sensitiveZ_[16] = 148.542;
	sensitiveZ_[17] = 158.536;
	sensitiveZ_[18] = 180.236;
	sensitiveZ_[19] = 190.23;
	sensitiveZ_[20] = 211.93;
	sensitiveZ_[21] = 221.924;
	sensitiveZ_[22] = 243.625;
	sensitiveZ_[23] = 253.619;
	sensitiveZ_[24] = 275.319;
	sensitiveZ_[25] = 285.313;
      }

      etaBoundary_.resize(indices_[3],0);
      etaBoundary_[0] =1.49188;
      etaBoundary_[1] =1.49329;
      etaBoundary_[2] =1.49475;
      etaBoundary_[3] =1.49614;
      etaBoundary_[4] =1.49758;
      etaBoundary_[5] =1.49894;
      etaBoundary_[6] =1.50036;
      etaBoundary_[7] =1.50171;
      etaBoundary_[8] =1.5031;
      etaBoundary_[9] =1.50443;
      etaBoundary_[10] =1.5058;
      etaBoundary_[11] =1.50711;
      etaBoundary_[12] =1.50846;
      etaBoundary_[13] =1.50974;
      etaBoundary_[14] =1.51108;
      etaBoundary_[15] =1.51234;
      etaBoundary_[16] =1.51365;
      etaBoundary_[17] =1.5149;
      etaBoundary_[18] =1.51619;
      etaBoundary_[19] =1.51743;
      etaBoundary_[20] =1.5187;
      etaBoundary_[21] =1.51991;
      etaBoundary_[22] =1.52116;
      etaBoundary_[23] =1.52236;
      etaBoundary_[24] =1.52359;
      etaBoundary_[25] =1.52821;
    } else if (versionNumber==83) {
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 26;
      indices_[4] = 48;
      indices_[5] = 52;
      indices_[6] = 62;

      sensitiveZ_.resize(indices_[6],0);
      sensitiveZ_[0] = 3213.95;
      sensitiveZ_[1] = 3223.94;
      sensitiveZ_[2] = 3245.64;
      sensitiveZ_[3] = 3255.64;
      sensitiveZ_[4] = 3277.34;
      sensitiveZ_[5] = 3287.33;
      sensitiveZ_[6] = 3309.03;
      sensitiveZ_[7] = 3319.03;
      sensitiveZ_[8] = 3340.73;
      sensitiveZ_[9] = 3350.72;
      sensitiveZ_[10] = 3372.42;
      sensitiveZ_[11] = 3382.42;
      sensitiveZ_[12] = 3404.12;
      sensitiveZ_[13] = 3414.11;
      sensitiveZ_[14] = 3435.81;
      sensitiveZ_[15] = 3445.8;
      sensitiveZ_[16] = 3467.5;
      sensitiveZ_[17] = 3477.5;
      sensitiveZ_[18] = 3499.2;
      sensitiveZ_[19] = 3509.19;
      sensitiveZ_[20] = 3530.89;
      sensitiveZ_[21] = 3540.89;
      sensitiveZ_[22] = 3562.59;
      sensitiveZ_[23] = 3572.58;
      sensitiveZ_[24] = 3594.28;
      sensitiveZ_[25] = 3604.28;
      sensitiveZ_[26] = 3671.83;
      sensitiveZ_[27] = 3734.88;
      sensitiveZ_[28] = 3797.93;
      sensitiveZ_[29] = 3860.98;
      sensitiveZ_[30] = 3924.03;
      sensitiveZ_[31] = 3987.08;
      sensitiveZ_[32] = 4050.13;
      sensitiveZ_[33] = 4113.18;
      sensitiveZ_[34] = 4176.23;
      sensitiveZ_[35] = 4239.28;
      sensitiveZ_[36] = 4302.33;
      sensitiveZ_[37] = 4384.58;
      sensitiveZ_[38] = 4466.83;
      sensitiveZ_[39] = 4549.08;
      sensitiveZ_[40] = 4631.33;
      sensitiveZ_[41] = 4713.58;
      sensitiveZ_[42] = 4795.83;
      sensitiveZ_[43] = 4878.08;
      sensitiveZ_[44] = 4960.33;
      sensitiveZ_[45] = 5042.58;
      sensitiveZ_[46] = 5124.83;
      sensitiveZ_[47] = 0.;
      sensitiveZ_[48] = 4109.88;
      sensitiveZ_[49] = 4172.93;
      sensitiveZ_[50] = 4235.98;
      sensitiveZ_[51] = 4299.03;
      sensitiveZ_[52] = 4381.28;
      sensitiveZ_[53] = 4463.53;
      sensitiveZ_[54] = 4545.78;
      sensitiveZ_[55] = 4628.03;
      sensitiveZ_[56] = 4710.28;
      sensitiveZ_[57] = 4792.53;
      sensitiveZ_[58] = 4874.78;
      sensitiveZ_[59] = 4957.03;
      sensitiveZ_[60] = 5039.28;
      sensitiveZ_[61] = 5121.53;

      etaBoundary_.resize(indices_[6],0);      
      etaBoundary_[0] =1.49188;
      etaBoundary_[1] =1.49329;
      etaBoundary_[2] =1.49475;
      etaBoundary_[3] =1.49614;
      etaBoundary_[4] =1.49758;
      etaBoundary_[5] =1.49894;
      etaBoundary_[6] =1.50036;
      etaBoundary_[7] =1.50171;
      etaBoundary_[8] =1.5031;
      etaBoundary_[9] =1.50443;
      etaBoundary_[10] =1.5058;
      etaBoundary_[11] =1.50711;
      etaBoundary_[12] =1.50846;
      etaBoundary_[13] =1.50974;
      etaBoundary_[14] =1.51108;
      etaBoundary_[15] =1.51234;
      etaBoundary_[16] =1.51365;
      etaBoundary_[17] =1.5149;
      etaBoundary_[18] =1.51619;
      etaBoundary_[19] =1.51743;
      etaBoundary_[20] =1.5187;
      etaBoundary_[21] =1.51991;
      etaBoundary_[22] =1.52116;
      etaBoundary_[23] =1.52236;
      etaBoundary_[24] =1.52359;
      etaBoundary_[25] =1.52477;
      etaBoundary_[26] =1.5298;
      etaBoundary_[27] =1.53436;
      etaBoundary_[28] =1.53879;
      etaBoundary_[29] =1.5431;
      etaBoundary_[30] =1.51781;
      etaBoundary_[31] =1.49332;
      etaBoundary_[32] =1.47031;
      etaBoundary_[33] =1.71063;
      etaBoundary_[34] =1.72487;
      etaBoundary_[35] =1.73892;
      etaBoundary_[36] =1.75279;
      etaBoundary_[37] =1.77061;
      etaBoundary_[38] =1.8869;
      etaBoundary_[39] =2.01578;
      etaBoundary_[40] =2.03305;
      etaBoundary_[41] =2.05005;
      etaBoundary_[42] =2.06676;
      etaBoundary_[43] =2.22908;
      etaBoundary_[44] =2.2454;
      etaBoundary_[45] =2.26146;
      etaBoundary_[46] =2.27727;
      etaBoundary_[47] =2.27727;
      etaBoundary_[48] =1.71063;
      etaBoundary_[49] =1.72487;
      etaBoundary_[50] =1.73892;
      etaBoundary_[51] =1.75279;
      etaBoundary_[52] =1.77061;
      etaBoundary_[53] =1.8869;
      etaBoundary_[54] =2.01578;
      etaBoundary_[55] =2.03305;
      etaBoundary_[56] =2.05005;
      etaBoundary_[57] =2.06676;
      etaBoundary_[58] =2.22908;
      etaBoundary_[59] =2.2454;
      etaBoundary_[60] =2.26146;
      etaBoundary_[61] =2.27727;
    }
    else if (versionNumber == 110){
      indices_[0] = 0;
      indices_[1] = 4;
      indices_[2] = indices_[1];
      indices_[3] = indices_[1];
      indices_[4] = indices_[1];
      indices_[5] = indices_[1];
      indices_[6] = indices_[1];
    }
    else {
      indices_[0] = 0;
      indices_[1] = versionNumber==24?11:10;
      indices_[2] = versionNumber==24?21:20;
      indices_[3] = versionNumber==24?31:30;
      indices_[4] = versionNumber==24?55:42;
      indices_[5] = versionNumber==24?65:54;
      indices_[6] = versionNumber==24?65:54;
    }
    
  };

  void buildDetector(const unsigned versionNumber,
		     const unsigned model=2,
		     bool concept=true,
		     bool isCaliceHcal=false,
		     bool bypassR=false);

  const HGCSSSubDetector & subDetectorByLayer(const unsigned aLayer);

  unsigned getSection(const unsigned aLayer) const;
  inline unsigned section(const DetectorEnum adet){
    if (enumMap_.find(adet) != enumMap_.end())
      return enumMap_[adet];
    return nSections();
  };

  void addSubdetector(const HGCSSSubDetector & adet);
  
  void finishInitialisation();

  inline bool isMixedLayer(const unsigned versionNumber,const unsigned aLayer){
    if (versionNumber!=63) return false;
    if (aLayer<36) return false;
    else return true;
  };

  inline unsigned nLayers(const unsigned aSection) const{
    return subdets_[aSection].nLayers();
  };

  inline unsigned nLayers(DetectorEnum adet){
    return subdets_[enumMap_[adet]].nLayers();
  };

  const HGCSSSubDetector & subDetectorByEnum(DetectorEnum adet);
  inline const HGCSSSubDetector & subDetectorBySection(const unsigned aSection) const{
    return subdets_[aSection];
  };

  inline unsigned nLayers() const{
    return nLayers_;
  };

  inline unsigned nSections() const{
    return nSections_;
  };

  inline DetectorEnum detType(const unsigned aSection) const{
    return subdets_[aSection].type;
  };
  
  inline DetectorEnum detTypeLayer(const unsigned aLayer) const{
    return subdets_[getSection(aLayer)].type;
  };
  
  inline std::string detName(const unsigned aSection) const{
    return subdets_[aSection].name;
  };


  void reset();

  void printDetector(std::ostream & aOs) const ;

  inline double sensitiveZ(const unsigned layer){
    if (layer<sensitiveZ_.size()) return sensitiveZ_[layer];
    else {
      std::cout << " ERROR! Trying to access layer " << layer << " outside of range: " << sensitiveZ_.size() << " nLayers " << nLayers_ << std::endl;
      exit(1);
    }
  };

  inline double etaBoundary(const unsigned layer){
    if (etaBoundary_.size()==0) return 0;
    if (layer<etaBoundary_.size()) return etaBoundary_[layer];
    else {
      std::cout << " ERROR! Trying to access layer " << layer << " outside of eta range: " << etaBoundary_.size() << " nLayers " << nLayers_ << std::endl;
      exit(1);
    }
  };

private:
  HGCSSDetector(){
    bypassRadius_ = false;
  };

  virtual ~HGCSSDetector(){
    reset();
  };
  
  std::vector<HGCSSSubDetector> subdets_;
  std::vector<unsigned> indices_;
  std::vector<unsigned> section_;
  std::map<DetectorEnum,unsigned> enumMap_;

  unsigned nLayers_;
  unsigned nSections_;
  bool bypassRadius_;

  std::vector<double> sensitiveZ_;
  std::vector<double> etaBoundary_;

};

HGCSSDetector & theDetector();




#endif
