/**
	commands 
**/
#define	INIT_DRAW 0
#define	REDRAW 1
#define	PSET 2
#define	PRELINE 3
#define	LINE 4
#define	CLS 5
#define	LINEATT 6
#define	ARCSTRUCT 7
#define	ARC 8
#define	INIT_SETP 10
#define	SETPOINTS 11
#define	PSETS 12
#define	FILL 13
#define	LINES 14
#define	BITMAP 15
#define	FILLOpaque 16
#define	CLIPMASK 17
#define	CLIPFULL 18
#define	ARROW 19
#define	SETARROW 20
#define	BITMAPS 21
#define	SETBITMAP 22
#define	WINDOW 23
#define	INIT_DRAW_RV 25
#define SETFONT 26
#define STRING 27
#define SETSTRING 28
#define STRINGS 29
#define FRAME 30
#define DIRECT 31
#define CYLINDER 32
#define CYLPSET 33
#define CYLLINE 34
#define CYLFRAME 35
#define	CLOSE 99
#define	WAIT 999
#define KEYIN 9996
#define STATICKEYIN 9995
#define	POINTER 9999
#define	CYLPOINT 9997
#define	POINT 9998
#define	MAGNIFY 1203
#define BITMAPOUT 500
#define MKARROW 501
/*
	set xhplot.init
*/
#define	AXES 1
#define	GRID 2
#define	ABNORMAL 3
#define	AXES_U 4
#define	AXES_R 5
#define	AXES_UR 6
#define	AXES_C 7
#define	OKDRAW 1
#define	NOTDRAW 2
#define	SETAGAIN 3
#define	NOLINE 4
#define	NOMOJI 5
/*
	line style
*/
#define	SOLID 1
#define	DASH 2
#define	DOTTED 3
#define	DOT_DASHED 4
#define	SHORT_DASHED 5
#define	LONG_DASHED 6
#define	ODD_DASHED 7
/*
	colors avairable
*/
#define	BLACK 0
#define	BLUE 1
#define	RED 2
#define	MAGENTA 3
#define	GREEN 4
#define	CYAN 5
#define	YELLOW 6
#define	WHITE 7
#define	STEELBLUE 8
#define	NAVYBLUE 9
#define	NAVY 10
#define	INDIANRED 11
#define	PINK 12
#define	DARKGREEN 13
#define	SPRINGGREEN 14
#define	TURQUOISE 15
/*
	arrow
*/
#define	HORIZONTAL 0
#define	VERTICAL 1
#define	TAIL 0
#define	HEAD 1
#define	NORMAL 0
#define	INVERSE 1

struct xddp_dat {
	int line_wid,line_att;
	int sv[4];
	double sc[4];
    double la[2];
    int lav[2];
    int fi[2];
    int kse;
    char font[100];
    int arrowtail[4];
    int flag;
} ;

struct xddp_dat xddp;

struct xddp_return {
	double x,y;
	int key;
};

struct xddp_return xddpret;

int xhplot(int mode,double x, double y, int icol);

void scalein(int *x, int *y);
void waitseconds(void);
void pointer(double xfactor, double yfactor,
				int sv[], double sc[], double *x, double *y);
void point(double xfactor, double yfactor,
				int sv[], double sc[], double *x, double *y);
void cylpoint(double xfactor, double yfactor,
				int sv[], double sc[], double *x, double *y);
void cliprectangle(int cp[]);
void setpoints(int ix, int iy, int n);
void clipmask(int n);
void linkpoints(int *n);
void fill(int flag, int n, int icol);
void setcolor(int icol);
void putpixel(int ix, int iy);
void putpixels(int n);
void setlineatt(void);
void line(int ix0,int iy0,int ix,int iy);
void lines(int n);
void setarcstruct(unsigned int *width, unsigned int *height,
					int *angle1, int *angle2);
void arc(int ix, int iy, unsigned int width,
			unsigned int height, int angle1,int angle2);
void putbitmap(void);
void gfcvt(double val, int ndigit, char s[]);
void moji(int ise, int jfx, int lbl, int lt, double t, int flag_ur);
void set_arrow(int *ver, int *head, int *inv);
void mkarrow(void);
void arrow(int ix, int iy, int ver, int head, int inv);
void draw_arrow(int i, int s, int ix, int iy);
int coordinate(int *ix, int *iy, double xfactor, double yfactor,
						double x, double y, int sv[], double sc[]);
int cylcoordinate(int *ix, int *iy, double xfactor, double yfactor,
						double x, double y, int sv[], double sc[]);
void setfont(void);
void drawstring(void);
void drawstrings(int x, int y, char string[]);
void setstring(char string[]);
void set_misc(int sv[], double sc[], double *xfactor,
				double *yfactor, double la[], int lav[],
					int arrowtail[]);
void struct2data(struct xddp_dat xddp, 
				int sv[], double sc[], double la[],
					int lav[], int fi[], int *kse, char *fontname,
						int arrowtail[], int *flag);

int no_wait_key(void);
int key_in(void);

void init(void);
double magx(int x);
double magy(int y);
double arcmagx(int x);
double arcmagy(int y);
double sgn(double x);
double arccoordinate(int ix, int iy, double xfactor, double yfactor,
			double *x, double *y, int sv[], double sc[]);
double arccylcoordinate(int ix, int iy, double xfactor, double yfactor,
			double *x, double *y, int sv[], double sc[]);

int isgn(double x);
int fiest(double x);

void data(int sv[], double sc[], double *xfactor,
		double *yfactor, double la[], int lav[], int arrowtail[]);
void axes_data(double la[],int lav[], int fi[], int *kse,
		int arrowtail[], int *flag);
void axes(int sv[], int lav[], int fi[], int *kse,
		double *xfactor, double *yfactor,
			double sc[], double la[], int icol, int flag);
void cylaxes(int sv[], int lav[], int fi[], int *kse,
		double *xfactor, double *yfactor, double sc[], double la[],
			int icol, int flag);

//int bitmaps(int xdest,int ydest,unsigned int width,
//		unsigned int height, Pixmap bitmap, int xhot, int yhot);
