#ifndef MIN
#define MIN(a, b) (a < b) ? a : b
#endif

#ifndef MAX
#define MAX(a, b) (a > b) ? a : b
#endif

#ifndef HSVtoRGB_h
#define HSVtoRGB_h

void RGBtoHSV( float r, float g, float b, float *h, float *s, float *v );
void HSVtoRGB( float h, float s, float v, float *r, float *g, float *b );

#endif
