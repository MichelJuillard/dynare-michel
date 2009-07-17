
#ifndef QT_H
#define QT_H

#define C_CHAR const char*
#define BLINT  int*
#define C_BLINT const int*
#define C_BLDOU const double*
#define BLDOU double*


extern "C"{
void ldm_(BLDOU, C_BLDOU, C_BLDOU, C_BLINT);
void ldsld_(BLDOU, C_BLDOU, C_BLDOU, C_BLINT);
void ldv_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
void mtt_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
void qt2ld_(BLDOU,C_BLDOU, C_BLINT);
void qt2t_(BLDOU, C_BLDOU, C_BLINT);
void s2d_(BLDOU,C_BLDOU, C_BLINT);
void s2u_(BLDOU,C_BLDOU, C_BLINT);
void td_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
void tm_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
void tstt_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
void tt_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
void tu_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
void tut_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
void tv_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
void qtv_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
#ifdef WINDOWS
void qtv_1__(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
#else
void qtv_1_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);
#endif
void qtsqtt_(BLDOU,C_BLDOU, C_BLDOU, C_BLINT);

};
#endif
