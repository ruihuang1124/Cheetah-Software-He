/* This file was automatically generated by CasADi 3.6.5.
 *  It consists of: 
 *   1) content generated by CasADi runtime: not copyrighted
 *   2) template code copied from CasADi source: permissively licensed (MIT-0)
 *   3) user code: owned by the user
 *
 */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) hkinodyn_casadi_ ## ID
#endif

#include <math.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

static const casadi_int casadi_s0[28] = {24, 1, 0, 24, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
static const casadi_int casadi_s1[5] = {1, 1, 0, 1, 0};
static const casadi_int casadi_s2[8] = {4, 1, 0, 4, 0, 1, 2, 3};

/* hkinodyn:(i0[24],i1[24],i2,i3[4])->(o0[24]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a4, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a5, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a6, a60, a61, a62, a63, a64, a65, a66, a67, a68, a69, a7, a70, a71, a72, a73, a74, a75, a76, a77, a78, a79, a8, a80, a9;
  a0=arg[0]? arg[0][0] : 0;
  a1=arg[0]? arg[0][2] : 0;
  a2=sin(a1);
  a3=arg[0]? arg[0][1] : 0;
  a4=cos(a3);
  a2=(a2/a4);
  a5=arg[0]? arg[0][7] : 0;
  a2=(a2*a5);
  a6=cos(a1);
  a6=(a6/a4);
  a7=arg[0]? arg[0][8] : 0;
  a6=(a6*a7);
  a2=(a2+a6);
  a6=arg[2]? arg[2][0] : 0;
  a2=(a2*a6);
  a2=(a0+a2);
  if (res[0]!=0) res[0][0]=a2;
  a2=cos(a1);
  a2=(a2*a5);
  a8=sin(a1);
  a8=(a8*a7);
  a2=(a2-a8);
  a2=(a2*a6);
  a2=(a3+a2);
  if (res[0]!=0) res[0][1]=a2;
  a2=arg[0]? arg[0][6] : 0;
  a8=sin(a1);
  a9=sin(a3);
  a8=(a8*a9);
  a8=(a8/a4);
  a8=(a8*a5);
  a8=(a2+a8);
  a9=cos(a1);
  a10=sin(a3);
  a9=(a9*a10);
  a9=(a9/a4);
  a9=(a9*a7);
  a8=(a8+a9);
  a8=(a8*a6);
  a8=(a1+a8);
  if (res[0]!=0) res[0][2]=a8;
  a8=arg[0]? arg[0][3] : 0;
  a9=arg[0]? arg[0][9] : 0;
  a4=(a9*a6);
  a4=(a8+a4);
  if (res[0]!=0) res[0][3]=a4;
  a4=arg[0]? arg[0][4] : 0;
  a10=arg[0]? arg[0][10] : 0;
  a11=(a10*a6);
  a11=(a4+a11);
  if (res[0]!=0) res[0][4]=a11;
  a11=arg[0]? arg[0][5] : 0;
  a12=arg[0]? arg[0][11] : 0;
  a13=(a12*a6);
  a13=(a11+a13);
  if (res[0]!=0) res[0][5]=a13;
  a13=1.1507568009567539e+01;
  a14=3.7947076036992650e-19;
  a15=(a14*a7);
  a16=-5.0200000000001471e-05;
  a17=(a16*a5);
  a15=(a15-a17);
  a15=(a15*a2);
  a17=2.6463700700179532e-01;
  a18=(a17*a7);
  a19=4.0822060971931784e-20;
  a20=(a19*a5);
  a18=(a18-a20);
  a18=(a18*a5);
  a15=(a15+a18);
  a18=(a19*a7);
  a20=2.9121873600000003e-01;
  a21=(a20*a5);
  a18=(a18-a21);
  a18=(a18*a7);
  a15=(a15+a18);
  a18=arg[3]? arg[3][0] : 0;
  a21=cos(a0);
  a22=sin(a3);
  a23=(a21*a22);
  a24=sin(a1);
  a25=(a23*a24);
  a0=sin(a0);
  a1=cos(a1);
  a26=(a0*a1);
  a25=(a25-a26);
  a26=arg[0]? arg[0][12] : 0;
  a27=(a26-a8);
  a28=(a25*a27);
  a29=(a21*a1);
  a30=(a0*a22);
  a31=(a30*a24);
  a29=(a29+a31);
  a31=arg[0]? arg[0][13] : 0;
  a32=(a31-a4);
  a33=(a29*a32);
  a28=(a28+a33);
  a3=cos(a3);
  a33=(a3*a24);
  a34=arg[0]? arg[0][14] : 0;
  a35=(a34-a11);
  a36=(a33*a35);
  a28=(a28+a36);
  a36=(a18*a28);
  a37=(a0*a24);
  a23=(a23*a1);
  a37=(a37+a23);
  a23=arg[1]? arg[1][0] : 0;
  a38=(a37*a23);
  a30=(a30*a1);
  a24=(a21*a24);
  a30=(a30-a24);
  a24=arg[1]? arg[1][1] : 0;
  a39=(a30*a24);
  a38=(a38+a39);
  a1=(a3*a1);
  a39=arg[1]? arg[1][2] : 0;
  a40=(a1*a39);
  a38=(a38+a40);
  a36=(a36*a38);
  a40=(a37*a27);
  a41=(a30*a32);
  a40=(a40+a41);
  a41=(a1*a35);
  a40=(a40+a41);
  a41=(a18*a40);
  a42=(a25*a23);
  a43=(a29*a24);
  a42=(a42+a43);
  a43=(a33*a39);
  a42=(a42+a43);
  a41=(a41*a42);
  a36=(a36-a41);
  a41=arg[3]? arg[3][1] : 0;
  a43=arg[0]? arg[0][15] : 0;
  a44=(a43-a8);
  a45=(a25*a44);
  a46=arg[0]? arg[0][16] : 0;
  a47=(a46-a4);
  a48=(a29*a47);
  a45=(a45+a48);
  a48=arg[0]? arg[0][17] : 0;
  a49=(a48-a11);
  a50=(a33*a49);
  a45=(a45+a50);
  a50=(a41*a45);
  a51=arg[1]? arg[1][3] : 0;
  a52=(a37*a51);
  a53=arg[1]? arg[1][4] : 0;
  a54=(a30*a53);
  a52=(a52+a54);
  a54=arg[1]? arg[1][5] : 0;
  a55=(a1*a54);
  a52=(a52+a55);
  a50=(a50*a52);
  a55=(a37*a44);
  a56=(a30*a47);
  a55=(a55+a56);
  a56=(a1*a49);
  a55=(a55+a56);
  a56=(a41*a55);
  a57=(a25*a51);
  a58=(a29*a53);
  a57=(a57+a58);
  a58=(a33*a54);
  a57=(a57+a58);
  a56=(a56*a57);
  a50=(a50-a56);
  a36=(a36+a50);
  a50=arg[3]? arg[3][2] : 0;
  a56=arg[0]? arg[0][18] : 0;
  a58=(a56-a8);
  a59=(a25*a58);
  a60=arg[0]? arg[0][19] : 0;
  a61=(a60-a4);
  a62=(a29*a61);
  a59=(a59+a62);
  a62=arg[0]? arg[0][20] : 0;
  a63=(a62-a11);
  a64=(a33*a63);
  a59=(a59+a64);
  a64=(a50*a59);
  a65=arg[1]? arg[1][6] : 0;
  a66=(a37*a65);
  a67=arg[1]? arg[1][7] : 0;
  a68=(a30*a67);
  a66=(a66+a68);
  a68=arg[1]? arg[1][8] : 0;
  a69=(a1*a68);
  a66=(a66+a69);
  a64=(a64*a66);
  a69=(a37*a58);
  a70=(a30*a61);
  a69=(a69+a70);
  a70=(a1*a63);
  a69=(a69+a70);
  a70=(a50*a69);
  a71=(a25*a65);
  a72=(a29*a67);
  a71=(a71+a72);
  a72=(a33*a68);
  a71=(a71+a72);
  a70=(a70*a71);
  a64=(a64-a70);
  a36=(a36+a64);
  a64=arg[3]? arg[3][3] : 0;
  a70=arg[0]? arg[0][21] : 0;
  a8=(a70-a8);
  a72=(a25*a8);
  a73=arg[0]? arg[0][22] : 0;
  a4=(a73-a4);
  a74=(a29*a4);
  a72=(a72+a74);
  a74=arg[0]? arg[0][23] : 0;
  a11=(a74-a11);
  a75=(a33*a11);
  a72=(a72+a75);
  a75=(a64*a72);
  a76=arg[1]? arg[1][9] : 0;
  a77=(a37*a76);
  a78=arg[1]? arg[1][10] : 0;
  a79=(a30*a78);
  a77=(a77+a79);
  a79=arg[1]? arg[1][11] : 0;
  a80=(a1*a79);
  a77=(a77+a80);
  a75=(a75*a77);
  a37=(a37*a8);
  a30=(a30*a4);
  a37=(a37+a30);
  a1=(a1*a11);
  a37=(a37+a1);
  a1=(a64*a37);
  a25=(a25*a76);
  a29=(a29*a78);
  a25=(a25+a29);
  a33=(a33*a79);
  a25=(a25+a33);
  a1=(a1*a25);
  a75=(a75-a1);
  a36=(a36+a75);
  a15=(a15+a36);
  a13=(a13*a15);
  a36=-1.6501345028411385e-17;
  a75=(a16*a2);
  a1=8.6899343001795346e-02;
  a33=(a1*a7);
  a75=(a75-a33);
  a75=(a75*a2);
  a33=(a19*a2);
  a29=(a14*a7);
  a33=(a33-a29);
  a33=(a33*a5);
  a75=(a75+a33);
  a20=(a20*a2);
  a33=(a16*a7);
  a20=(a20-a33);
  a20=(a20*a7);
  a75=(a75+a20);
  a40=(a18*a40);
  a21=(a21*a3);
  a20=(a21*a23);
  a0=(a0*a3);
  a3=(a0*a24);
  a20=(a20+a3);
  a3=(a22*a39);
  a20=(a20-a3);
  a40=(a40*a20);
  a27=(a21*a27);
  a32=(a0*a32);
  a27=(a27+a32);
  a35=(a22*a35);
  a27=(a27-a35);
  a35=(a18*a27);
  a35=(a35*a38);
  a40=(a40-a35);
  a55=(a41*a55);
  a35=(a21*a51);
  a38=(a0*a53);
  a35=(a35+a38);
  a38=(a22*a54);
  a35=(a35-a38);
  a55=(a55*a35);
  a44=(a21*a44);
  a47=(a0*a47);
  a44=(a44+a47);
  a49=(a22*a49);
  a44=(a44-a49);
  a49=(a41*a44);
  a49=(a49*a52);
  a55=(a55-a49);
  a40=(a40+a55);
  a69=(a50*a69);
  a55=(a21*a65);
  a49=(a0*a67);
  a55=(a55+a49);
  a49=(a22*a68);
  a55=(a55-a49);
  a69=(a69*a55);
  a58=(a21*a58);
  a61=(a0*a61);
  a58=(a58+a61);
  a63=(a22*a63);
  a58=(a58-a63);
  a63=(a50*a58);
  a63=(a63*a66);
  a69=(a69-a63);
  a40=(a40+a69);
  a37=(a64*a37);
  a69=(a21*a76);
  a63=(a0*a78);
  a69=(a69+a63);
  a63=(a22*a79);
  a69=(a69-a63);
  a37=(a37*a69);
  a21=(a21*a8);
  a0=(a0*a4);
  a21=(a21+a0);
  a22=(a22*a11);
  a21=(a21-a22);
  a22=(a64*a21);
  a22=(a22*a77);
  a37=(a37-a22);
  a40=(a40+a37);
  a75=(a75+a40);
  a40=(a36*a75);
  a13=(a13+a40);
  a40=1.9836632835337469e-03;
  a1=(a1*a5);
  a37=(a14*a2);
  a1=(a1-a37);
  a1=(a1*a2);
  a14=(a14*a5);
  a17=(a17*a2);
  a14=(a14-a17);
  a14=(a14*a5);
  a1=(a1+a14);
  a16=(a16*a5);
  a19=(a19*a2);
  a16=(a16-a19);
  a16=(a16*a7);
  a1=(a1+a16);
  a27=(a18*a27);
  a27=(a27*a42);
  a28=(a18*a28);
  a28=(a28*a20);
  a27=(a27-a28);
  a44=(a41*a44);
  a44=(a44*a57);
  a45=(a41*a45);
  a45=(a45*a35);
  a44=(a44-a45);
  a27=(a27+a44);
  a58=(a50*a58);
  a58=(a58*a71);
  a59=(a50*a59);
  a59=(a59*a55);
  a58=(a58-a59);
  a27=(a27+a58);
  a21=(a64*a21);
  a21=(a21*a25);
  a72=(a64*a72);
  a72=(a72*a69);
  a21=(a21-a72);
  a27=(a27+a21);
  a1=(a1+a27);
  a27=(a40*a1);
  a13=(a13+a27);
  a13=(a13*a6);
  a2=(a2+a13);
  if (res[0]!=0) res[0][6]=a2;
  a36=(a36*a15);
  a2=3.7787609954083861e+00;
  a2=(a2*a75);
  a36=(a36+a2);
  a2=-5.3253846714500212e-19;
  a13=(a2*a1);
  a36=(a36+a13);
  a36=(a36*a6);
  a5=(a5+a36);
  if (res[0]!=0) res[0][7]=a5;
  a40=(a40*a15);
  a2=(a2*a75);
  a40=(a40+a2);
  a2=3.4338453401566058e+00;
  a2=(a2*a1);
  a40=(a40+a2);
  a40=(a40*a6);
  a7=(a7+a40);
  if (res[0]!=0) res[0][8]=a7;
  a23=(a18*a23);
  a51=(a41*a51);
  a23=(a23+a51);
  a65=(a50*a65);
  a23=(a23+a65);
  a76=(a64*a76);
  a23=(a23+a76);
  a76=8.9119999999999990e+00;
  a23=(a23/a76);
  a23=(a23*a6);
  a9=(a9+a23);
  if (res[0]!=0) res[0][9]=a9;
  a24=(a18*a24);
  a53=(a41*a53);
  a24=(a24+a53);
  a67=(a50*a67);
  a24=(a24+a67);
  a78=(a64*a78);
  a24=(a24+a78);
  a24=(a24/a76);
  a24=(a24*a6);
  a10=(a10+a24);
  if (res[0]!=0) res[0][10]=a10;
  a10=-9.8100000000000005e+00;
  a39=(a18*a39);
  a54=(a41*a54);
  a39=(a39+a54);
  a68=(a50*a68);
  a39=(a39+a68);
  a79=(a64*a79);
  a39=(a39+a79);
  a39=(a39/a76);
  a10=(a10+a39);
  a10=(a10*a6);
  a12=(a12+a10);
  if (res[0]!=0) res[0][11]=a12;
  a12=5.0000000000000000e-01;
  a18=(a12<=a18);
  a18=(!a18);
  a10=arg[1]? arg[1][12] : 0;
  a10=(a10*a6);
  a10=(a18?a10:0);
  a26=(a26+a10);
  if (res[0]!=0) res[0][12]=a26;
  a26=arg[1]? arg[1][13] : 0;
  a26=(a26*a6);
  a26=(a18?a26:0);
  a31=(a31+a26);
  if (res[0]!=0) res[0][13]=a31;
  a31=arg[1]? arg[1][14] : 0;
  a31=(a31*a6);
  a18=(a18?a31:0);
  a34=(a34+a18);
  if (res[0]!=0) res[0][14]=a34;
  a41=(a12<=a41);
  a41=(!a41);
  a34=arg[1]? arg[1][15] : 0;
  a34=(a34*a6);
  a34=(a41?a34:0);
  a43=(a43+a34);
  if (res[0]!=0) res[0][15]=a43;
  a43=arg[1]? arg[1][16] : 0;
  a43=(a43*a6);
  a43=(a41?a43:0);
  a46=(a46+a43);
  if (res[0]!=0) res[0][16]=a46;
  a46=arg[1]? arg[1][17] : 0;
  a46=(a46*a6);
  a41=(a41?a46:0);
  a48=(a48+a41);
  if (res[0]!=0) res[0][17]=a48;
  a50=(a12<=a50);
  a50=(!a50);
  a48=arg[1]? arg[1][18] : 0;
  a48=(a48*a6);
  a48=(a50?a48:0);
  a56=(a56+a48);
  if (res[0]!=0) res[0][18]=a56;
  a56=arg[1]? arg[1][19] : 0;
  a56=(a56*a6);
  a56=(a50?a56:0);
  a60=(a60+a56);
  if (res[0]!=0) res[0][19]=a60;
  a60=arg[1]? arg[1][20] : 0;
  a60=(a60*a6);
  a50=(a50?a60:0);
  a62=(a62+a50);
  if (res[0]!=0) res[0][20]=a62;
  a12=(a12<=a64);
  a12=(!a12);
  a64=arg[1]? arg[1][21] : 0;
  a64=(a64*a6);
  a64=(a12?a64:0);
  a70=(a70+a64);
  if (res[0]!=0) res[0][21]=a70;
  a70=arg[1]? arg[1][22] : 0;
  a70=(a70*a6);
  a70=(a12?a70:0);
  a73=(a73+a70);
  if (res[0]!=0) res[0][22]=a73;
  a73=arg[1]? arg[1][23] : 0;
  a73=(a73*a6);
  a12=(a12?a73:0);
  a74=(a74+a12);
  if (res[0]!=0) res[0][23]=a74;
  return 0;
}

CASADI_SYMBOL_EXPORT int hkinodyn(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int hkinodyn_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int hkinodyn_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void hkinodyn_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int hkinodyn_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void hkinodyn_release(int mem) {
}

CASADI_SYMBOL_EXPORT void hkinodyn_incref(void) {
}

CASADI_SYMBOL_EXPORT void hkinodyn_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int hkinodyn_n_in(void) { return 4;}

CASADI_SYMBOL_EXPORT casadi_int hkinodyn_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real hkinodyn_default_in(casadi_int i) {
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* hkinodyn_name_in(casadi_int i) {
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* hkinodyn_name_out(casadi_int i) {
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* hkinodyn_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s1;
    case 3: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* hkinodyn_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int hkinodyn_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 4;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

CASADI_SYMBOL_EXPORT int hkinodyn_work_bytes(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 4*sizeof(const casadi_real*);
  if (sz_res) *sz_res = 1*sizeof(casadi_real*);
  if (sz_iw) *sz_iw = 0*sizeof(casadi_int);
  if (sz_w) *sz_w = 0*sizeof(casadi_real);
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
