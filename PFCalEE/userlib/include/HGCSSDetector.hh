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
  ~HGCSSSubDetector(){};

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
    else if (versionNumber == 30 || versionNumber == 60 || (versionNumber >= 100 && versionNumber < 104)){
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
	sensitiveZ_[0] = -54.3;
	sensitiveZ_[1] = -45.2;
	sensitiveZ_[2] = -29.9;
	sensitiveZ_[3] = -20.8;
	sensitiveZ_[4] = -5.5;
	sensitiveZ_[5] = 3.6;
	sensitiveZ_[6] = 18.9;
	sensitiveZ_[7] = 28;
	sensitiveZ_[8] = 43.3;
	sensitiveZ_[9] = 52.4;
	sensitiveZ_[10] = 67.7;
	sensitiveZ_[11] = 76.8;
	sensitiveZ_[12] = 92.1;
	sensitiveZ_[13] = 101.2;
	sensitiveZ_[14] = 116.5;
	sensitiveZ_[15] = 125.6;
	sensitiveZ_[16] = 140.9;
	sensitiveZ_[17] = 150;
	sensitiveZ_[18] = 165.3;
	sensitiveZ_[19] = 174.4;
	sensitiveZ_[20] = 189.7;
	sensitiveZ_[21] = 198.8;
	sensitiveZ_[22] = 214.1;
	sensitiveZ_[23] = 223.2;
	sensitiveZ_[24] = 238.5;
	sensitiveZ_[25] = 247.6;
	sensitiveZ_[26] = 262.9;
	sensitiveZ_[27] = 272;
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
      sensitiveZ_[0] = -73.05;
      sensitiveZ_[1] = -63.95;
      sensitiveZ_[2] = -45.65;
      sensitiveZ_[3] = -36.55;
      sensitiveZ_[4] = -18.25;
      sensitiveZ_[5] = -9.15;
      sensitiveZ_[6] = 9.15;
      sensitiveZ_[7] = 18.25;
      sensitiveZ_[8] = 36.55;
      sensitiveZ_[9] = 45.65;
      sensitiveZ_[10] = 63.95;
      sensitiveZ_[11] = 73.05;
      sensitiveZ_[12] = 91.35;
      sensitiveZ_[13] = 100.45;
      sensitiveZ_[14] = 118.75;
      sensitiveZ_[15] = 127.85;
      sensitiveZ_[16] = 146.15;
      sensitiveZ_[17] = 155.25;
      sensitiveZ_[18] = 173.55;
      sensitiveZ_[19] = 182.65;
      sensitiveZ_[20] = 200.95;
      sensitiveZ_[21] = 210.05;
      sensitiveZ_[22] = 228.35;
      sensitiveZ_[23] = 237.45;
      sensitiveZ_[24] = 255.75;
      sensitiveZ_[25] = 264.85;
      sensitiveZ_[26] = 283.15;
      sensitiveZ_[27] = 292.25;
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
	sensitiveZ_[0] = -85.55;
	sensitiveZ_[1] = -76.45;
	sensitiveZ_[2] = -56.15;
	sensitiveZ_[3] = -47.05;
	sensitiveZ_[4] = -26.75;
	sensitiveZ_[5] = -17.65;
	sensitiveZ_[6] = 2.65;
	sensitiveZ_[7] = 11.75;
	sensitiveZ_[8] = 32.05;
	sensitiveZ_[9] = 41.15;
	sensitiveZ_[10] = 61.45;
	sensitiveZ_[11] = 70.55;
	sensitiveZ_[12] = 90.85;
	sensitiveZ_[13] = 99.95;
	sensitiveZ_[14] = 120.25;
	sensitiveZ_[15] = 129.35;
	sensitiveZ_[16] = 149.65;
	sensitiveZ_[17] = 158.75;
	sensitiveZ_[18] = 179.05;
	sensitiveZ_[19] = 188.15;
	sensitiveZ_[20] = 208.45;
	sensitiveZ_[21] = 217.55;
	sensitiveZ_[22] = 237.85;
	sensitiveZ_[23] = 246.95;
	sensitiveZ_[24] = 267.25;
	sensitiveZ_[25] = 276.35;
	sensitiveZ_[26] = 296.65;
	sensitiveZ_[27] = 305.75;
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
      //for (unsigned iL(0); iL<36; ++iL){
      //etaBoundary_[iL] = 1.4;
      //}
      // outer eta boundary for layers 0-35
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

  ~HGCSSDetector(){
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
