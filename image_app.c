// image_app.c
// Single-file PGM image processing application (supports P2 and P5)
// Includes bilinear zoom, shrink (block average), filters, edges, LBP, diagnostics
// Compile: gcc -O2 -std=c99 -lm image_app.c -o image_app

#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

typedef struct {
    int width;
    int height;
    int maxval;       // internally 255
    unsigned char *data; // length = width * height, grayscale
} Image;

/* ------------------ Utility ------------------ */
static void die(const char *msg) { fprintf(stderr, "Error: %s\n", msg); exit(EXIT_FAILURE); }
Image *create_image(int w, int h, int maxv) {
    Image *img = (Image*)malloc(sizeof(Image));
    if (!img) die("out of memory");
    img->width = w; img->height = h; img->maxval = maxv;
    img->data = (unsigned char*)calloc((size_t)w * h, sizeof(unsigned char));
    if (!img->data) die("out of memory");
    return img;
}
void free_image(Image *img) { if (!img) return; free(img->data); free(img); }

/* skip whitespace and comments, return next non-space char (stream left at it) */
static int next_nonspace(FILE *f) {
    int c;
    while ((c = fgetc(f)) != EOF) {
        if (isspace(c)) continue;
        if (c == '#') {
            while ((c = fgetc(f)) != EOF && c != '\n');
            continue;
        }
        ungetc(c, f);
        return c;
    }
    return EOF;
}

/* ------------------ PGM I/O ------------------ */
Image *load_pgm(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;
    char magic[3] = {0};
    if (fscanf(f, "%2s", magic) != 1) { fclose(f); return NULL; }
    if (strcmp(magic, "P2") != 0 && strcmp(magic, "P5") != 0) { fclose(f); return NULL; }
    int ascii = (strcmp(magic, "P2") == 0);
    if (next_nonspace(f) == EOF) { fclose(f); return NULL; }
    int w, h, maxv;
    if (fscanf(f, "%d", &w) != 1) { fclose(f); return NULL; }
    if (next_nonspace(f) == EOF) { fclose(f); return NULL; }
    if (fscanf(f, "%d", &h) != 1) { fclose(f); return NULL; }
    if (next_nonspace(f) == EOF) { fclose(f); return NULL; }
    if (fscanf(f, "%d", &maxv) != 1) { fclose(f); return NULL; }
    if (maxv <= 0 || maxv > 65535) { fclose(f); return NULL; }
    int c = fgetc(f);
    if (c == '\r') { int c2 = fgetc(f); if (c2 != '\n') ungetc(c2, f); }
    else if (c != '\n' && !isspace(c)) ungetc(c, f);
    Image *img = create_image(w, h, maxv);
    size_t npix = (size_t)w * h;
    if (ascii) {
        for (size_t i=0; i < npix; ++i) {
            int v; if (fscanf(f, "%d", &v) != 1) { free_image(img); fclose(f); return NULL; }
            if (v < 0) v = 0; if (v > img->maxval) v = img->maxval;
            if (img->maxval == 255) img->data[i] = (unsigned char)v;
            else img->data[i] = (unsigned char) round((double)v * 255.0 / img->maxval);
        }
    } else {
        if (img->maxval < 256) {
            size_t read = fread(img->data, 1, npix, f);
            if (read != npix) { free_image(img); fclose(f); return NULL; }
        } else {
            for (size_t i=0; i < npix; ++i) {
                int hi = fgetc(f); int lo = fgetc(f);
                if (hi == EOF || lo == EOF) { free_image(img); fclose(f); return NULL; }
                int v = (hi << 8) | lo;
                if (v < 0) v = 0; if (v > img->maxval) v = img->maxval;
                img->data[i] = (unsigned char) round((double)v * 255.0 / img->maxval);
            }
        }
    }
    fclose(f);
    img->maxval = 255;
    return img;
}

int save_pgm(const Image *img, const char *path, int ascii) {
    if (!img) return 0;
    FILE *f = fopen(path, "wb");
    if (!f) return 0;
    if (ascii) {
        fprintf(f, "P2\n# Created by image_app\n%d %d\n255\n", img->width, img->height);
        size_t np = (size_t)img->width * img->height;
        int perline = 17;
        for (size_t i = 0; i < np; ++i) {
            fprintf(f, "%d", img->data[i]);
            if ((i % perline) == perline-1 || i == np-1) fprintf(f, "\n");
            else fprintf(f, " ");
        }
    } else {
        fprintf(f, "P5\n# Created by image_app\n%d %d\n255\n", img->width, img->height);
        size_t np = (size_t)img->width * img->height;
        size_t written = fwrite(img->data, 1, np, f);
        if (written != np) { fclose(f); return 0; }
    }
    fclose(f);
    return 1;
}

/* ------------------ Basic ops ------------------ */
Image *duplicate_image(const Image *src) {
    Image *dst = create_image(src->width, src->height, src->maxval);
    memcpy(dst->data, src->data, (size_t)src->width * src->height);
    return dst;
}

/* Bilinear zoom by integer factor (>=1) */
Image *zoom_bilinear(const Image *src, int factor) {
    if (factor <= 0) return NULL;
    if (factor == 1) return duplicate_image(src);
    int w_out = src->width * factor;
    int h_out = src->height * factor;
    Image *dst = create_image(w_out, h_out, 255);
    for (int y_out = 0; y_out < h_out; ++y_out) {
        double gy = (double)y_out / (double)factor;
        int y0 = (int)floor(gy);
        int y1 = y0 + 1;
        double wy = gy - y0;
        if (y0 < 0) y0 = 0;
        if (y1 >= src->height) y1 = src->height - 1;
        for (int x_out = 0; x_out < w_out; ++x_out) {
            double gx = (double)x_out / (double)factor;
            int x0 = (int)floor(gx);
            int x1 = x0 + 1;
            double wx = gx - x0;
            if (x0 < 0) x0 = 0;
            if (x1 >= src->width) x1 = src->width - 1;
            // fetch four neighbors
            unsigned char I00 = src->data[y0*src->width + x0];
            unsigned char I10 = src->data[y0*src->width + x1];
            unsigned char I01 = src->data[y1*src->width + x0];
            unsigned char I11 = src->data[y1*src->width + x1];
            double val = (1-wx)*(1-wy)*I00 + wx*(1-wy)*I10 + (1-wx)*wy*I01 + wx*wy*I11;
            int iv = (int)round(val);
            if (iv < 0) iv = 0; if (iv > 255) iv = 255;
            dst->data[y_out*w_out + x_out] = (unsigned char)iv;
        }
    }
    return dst;
}

/* Shrink by integer factor using block average */
Image *shrink_avg(const Image *src, int factor) {
    if (factor <= 0) return NULL;
    if (factor == 1) return duplicate_image(src);
    int w = src->width / factor; if (w < 1) w = 1;
    int h = src->height / factor; if (h < 1) h = 1;
    Image *dst = create_image(w, h, 255);
    for (int y=0; y<h; ++y) {
        for (int x=0; x<w; ++x) {
            int sx = x * factor;
            int sy = y * factor;
            int sum = 0, cnt = 0;
            for (int dy=0; dy<factor; ++dy) {
                int yy = sy + dy; if (yy >= src->height) break;
                for (int dx=0; dx<factor; ++dx) {
                    int xx = sx + dx; if (xx >= src->width) break;
                    sum += src->data[yy*src->width + xx]; cnt++;
                }
            }
            dst->data[y*w + x] = (unsigned char) round((double)sum / (cnt?cnt:1));
        }
    }
    return dst;
}

/* ------------------ Filters (box, gaussian, median) ------------------ */

Image *box_filter(const Image *src, int ksize) {
    if (ksize <= 1) return duplicate_image(src);
    if (ksize % 2 == 0) ksize++;
    int r = ksize/2; int w = src->width, h = src->height;
    Image *dst = create_image(w,h,255);
    for (int y=0;y<h;y++){
        for (int x=0;x<w;x++){
            int sum=0, cnt=0;
            for (int dy=-r; dy<=r; ++dy){
                int yy=y+dy; if (yy<0||yy>=h) continue;
                for (int dx=-r; dx<=r; ++dx){
                    int xx=x+dx; if (xx<0||xx>=w) continue;
                    sum += src->data[yy*w + xx]; cnt++;
                }
            }
            dst->data[y*w + x] = (unsigned char) round((double)sum/(cnt?cnt:1));
        }
    }
    return dst;
}

Image *gaussian_filter(const Image *src, int ksize, double sigma) {
    if (ksize <= 1) return duplicate_image(src);
    if (ksize % 2 == 0) ksize++;
    if (sigma <= 0) sigma = 1.0;
    int r = ksize/2; double *kernel = (double*)malloc(sizeof(double)*ksize*ksize);
    double sumk = 0.0;
    for (int y=-r;y<=r;y++){
        for (int x=-r;x<=r;x++){
            double v = exp(-(x*x + y*y)/(2*sigma*sigma));
            kernel[(y+r)*ksize + (x+r)] = v; sumk += v;
        }
    }
    for (int i=0;i<ksize*ksize;i++) kernel[i]/=sumk;
    int w = src->width, h = src->height;
    Image *dst = create_image(w,h,255);
    for (int y=0;y<h;y++){
        for (int x=0;x<w;x++){
            double acc=0.0;
            for (int dy=-r;dy<=r;dy++){
                int yy=y+dy; if (yy<0||yy>=h) continue;
                for (int dx=-r;dx<=r;dx++){
                    int xx=x+dx; if (xx<0||xx>=w) continue;
                    acc += src->data[yy*w + xx] * kernel[(dy+r)*ksize + (dx+r)];
                }
            }
            int v = (int)round(acc); if (v<0)v=0;if(v>255)v=255;
            dst->data[y*w + x] = (unsigned char)v;
        }
    }
    free(kernel); return dst;
}

Image *median_filter(const Image *src, int ksize) {
    if (ksize <= 1) return duplicate_image(src);
    if (ksize % 2 == 0) ksize++;
    int r = ksize/2; int w = src->width, h = src->height;
    Image *dst = create_image(w,h,255);
    int maxcells = ksize*ksize; unsigned char *buf = (unsigned char*)malloc(maxcells);
    for (int y=0;y<h;y++){
        for (int x=0;x<w;x++){
            int cnt=0;
            for (int dy=-r;dy<=r;dy++){
                int yy=y+dy; if (yy<0||yy>=h) continue;
                for (int dx=-r;dx<=r;dx++){
                    int xx=x+dx; if (xx<0||xx>=w) continue;
                    buf[cnt++]=src->data[yy*w + xx];
                }
            }
            for (int i=1;i<cnt;i++){ unsigned char key=buf[i]; int j=i-1; while(j>=0&&buf[j]>key){buf[j+1]=buf[j]; j--; } buf[j+1]=key; }
            dst->data[y*w + x] = buf[cnt/2];
        }
    }
    free(buf); return dst;
}

/* ------------------ Edge detectors ------------------ */
/* (Sobel, Prewitt, Canny — same approach as earlier) */

static void sobel_gradients(const Image *src, double *gx, double *gy, double *gmag) {
    int w = src->width, h = src->height;
    int Gx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
    int Gy[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};
    for (int y=0;y<h;y++){
        for (int x=0;x<w;x++){
            double xacc=0.0, yacc=0.0;
            for (int dy=-1;dy<=1;dy++){
                int yy=y+dy; if (yy<0||yy>=h) continue;
                for (int dx=-1;dx<=1;dx++){
                    int xx=x+dx; if (xx<0||xx>=w) continue;
                    unsigned char p = src->data[yy*w + xx];
                    xacc += Gx[dy+1][dx+1] * p;
                    yacc += Gy[dy+1][dx+1] * p;
                }
            }
            int idx = y*w + x;
            gx[idx] = xacc; gy[idx] = yacc; gmag[idx] = hypot(xacc, yacc);
        }
    }
}

Image *sobel_edge(const Image *src) {
    int w = src->width, h = src->height;
    double *gx = (double*)calloc((size_t)w*h, sizeof(double));
    double *gy = (double*)calloc((size_t)w*h, sizeof(double));
    double *gmag = (double*)calloc((size_t)w*h, sizeof(double));
    sobel_gradients(src, gx, gy, gmag);
    double maxv = 0.0; for (int i=0;i<w*h;i++) if (gmag[i] > maxv) maxv = gmag[i]; if (maxv<=0) maxv=1.0;
    Image *dst = create_image(w,h,255);
    for (int i=0;i<w*h;i++) { int v = (int)round(255.0 * gmag[i] / maxv); if (v<0)v=0; if (v>255)v=255; dst->data[i] = (unsigned char)v; }
    free(gx); free(gy); free(gmag); return dst;
}

Image *prewitt_edge(const Image *src) {
    int w = src->width, h = src->height;
    double *gmag = (double*)calloc((size_t)w*h, sizeof(double));
    int Gx[3][3] = {{-1,0,1},{-1,0,1},{-1,0,1}};
    int Gy[3][3] = {{1,1,1},{0,0,0},{-1,-1,-1}};
    for (int y=0;y<h;y++){
        for (int x=0;x<w;x++){
            double xacc=0.0, yacc=0.0;
            for (int dy=-1;dy<=1;dy++){
                int yy=y+dy; if (yy<0||yy>=h) continue;
                for (int dx=-1;dx<=1;dx++){
                    int xx=x+dx; if (xx<0||xx>=w) continue;
                    unsigned char p = src->data[yy*w + xx];
                    xacc += Gx[dy+1][dx+1] * p; yacc += Gy[dy+1][dx+1] * p;
                }
            }
            gmag[y*w + x] = hypot(xacc, yacc);
        }
    }
    double maxv = 0.0; for (int i=0;i<w*h;i++) if (gmag[i] > maxv) maxv = gmag[i]; if (maxv<=0) maxv=1.0;
    Image *dst = create_image(w,h,255);
    for (int i=0;i<w*h;i++) { int v = (int)round(255.0 * gmag[i] / maxv); if (v<0)v=0; if (v>255)v=255; dst->data[i] = (unsigned char)v; }
    free(gmag); return dst;
}

Image *canny_edge(const Image *src, double low_ratio, double high_ratio) {
    if (low_ratio < 0) low_ratio = 0.05; if (high_ratio <= low_ratio) high_ratio = low_ratio * 3;
    Image *smooth = gaussian_filter(src, 5, 1.0);
    int w = smooth->width, h = smooth->height;
    double *gx = (double*)calloc((size_t)w*h, sizeof(double));
    double *gy = (double*)calloc((size_t)w*h, sizeof(double));
    double *gmag = (double*)calloc((size_t)w*h, sizeof(double));
    sobel_gradients(smooth, gx, gy, gmag);
    double *nms = (double*)calloc((size_t)w*h, sizeof(double));
    for (int y=1;y<h-1;y++){
        for (int x=1;x<w-1;x++){
            int idx = y*w + x;
            double angle = atan2(gy[idx], gx[idx]) * 180.0 / M_PI; if (angle < 0) angle += 180.0;
            double m = gmag[idx]; double m1=0,m2=0;
            if ((angle >=0 && angle <22.5) || (angle >=157.5 && angle <=180)) { m1=gmag[idx+1]; m2=gmag[idx-1]; }
            else if (angle >=22.5 && angle <67.5) { m1=gmag[idx - w + 1]; m2=gmag[idx + w -1]; }
            else if (angle >=67.5 && angle <112.5) { m1=gmag[idx - w]; m2=gmag[idx + w]; }
            else { m1=gmag[idx - w -1]; m2=gmag[idx + w +1]; }
            nms[idx] = (m >= m1 && m >= m2) ? m : 0;
        }
    }
    double maxv = 0.0; for (int i=0;i<w*h;i++) if (nms[i] > maxv) maxv = nms[i]; if (maxv<=0) maxv=1.0;
    double high = high_ratio * maxv; double low = low_ratio * maxv;
    unsigned char *mark = (unsigned char*)calloc((size_t)w*h, sizeof(unsigned char));
    enum { NONE=0, WEAK=1, STRONG=2 };
    for (int i=0;i<w*h;i++) { if (nms[i] >= high) mark[i]=STRONG; else if (nms[i] >= low) mark[i]=WEAK; else mark[i]=NONE; }
    int changed = 1; while (changed) { changed = 0; for (int y=1;y<h-1;y++){ for (int x=1;x<w-1;x++){ int idx=y*w+x; if (mark[idx]==WEAK){ int found=0; for (int dy=-1;dy<=1 && !found;dy++) for (int dx=-1;dx<=1;dx++){ if (dy==0 && dx==0) continue; int nidx=(y+dy)*w+(x+dx); if (mark[nidx]==STRONG){found=1; break;} } if (found) { mark[idx]=STRONG; changed=1; } } } } }
    Image *dst = create_image(w,h,255);
    for (int i=0;i<w*h;i++) dst->data[i] = (mark[i]==STRONG) ? 255 : 0;
    free_image(smooth); free(gx); free(gy); free(gmag); free(nms); free(mark);
    return dst;
}

/* ------------------ LBP ------------------ */
Image *compute_lbp(const Image *src) {
    int w = src->width, h = src->height;
    Image *dst = create_image(w,h,255);
    for (int y=1;y<h-1;y++){
        for (int x=1;x<w-1;x++){
            unsigned char center = src->data[y*w + x];
            unsigned char code = 0;
            code |= (src->data[(y-1)*w + (x-1)] >= center) << 7;
            code |= (src->data[(y-1)*w + (x  )] >= center) << 6;
            code |= (src->data[(y-1)*w + (x+1)] >= center) << 5;
            code |= (src->data[(y  )*w + (x+1)] >= center) << 4;
            code |= (src->data[(y+1)*w + (x+1)] >= center) << 3;
            code |= (src->data[(y+1)*w + (x  )] >= center) << 2;
            code |= (src->data[(y+1)*w + (x-1)] >= center) << 1;
            code |= (src->data[(y  )*w + (x-1)] >= center) << 0;
            dst->data[y*w + x] = code;
        }
    }
    return dst;
}

/* ------------------ Diagnostics ------------------ */
/* Returns sum of all pixel values (64-bit) */
unsigned long long image_checksum(const Image *img) {
    unsigned long long s = 0;
    size_t np = (size_t)img->width * img->height;
    for (size_t i=0;i<np;i++) s += img->data[i];
    return s;
}

/* Print top-left 5x5 pixel block (or smaller if image smaller) */
void print_sample(const Image *img) {
    int W = img->width, H = img->height;
    int sw = (W < 5) ? W : 5;
    int sh = (H < 5) ? H : 5;
    printf("Top-left %dx%d sample:\n", sw, sh);
    for (int y=0;y<sh;y++){
        for (int x=0;x<sw;x++){
            printf("%3d ", img->data[y*W + x]);
        }
        printf("\n");
    }
}

/* ------------------ Menu / Main ------------------ */
static void print_menu(void) {
    puts("\n===== Image Processing Application (PGM only) =====");
    puts("1 - Load a PGM image file (P2 or P5)");
    puts("2 - Zoom / Shrink the loaded image by integer factor");
    puts("3 - Apply Filter on the Image");
    puts("    3.1 Average (box)");
    puts("    3.2 Mean (Gaussian)");
    puts("    3.3 Median");
    puts("4 - Edge Detection");
    puts("    4.1 Canny");
    puts("    4.2 Sobel");
    puts("    4.3 Prewitt");
    puts("5 - Compute and store LBP code for each pixel");
    puts("6 - Save the processed image (choose P2 or P5)");
    puts("7 - Show current image info");
    puts("8 - Diagnostic: checksum + sample (helps verify zoom)");
    puts("0 - Exit");
    printf("Enter option: ");
}

static void read_line_discard(void) { int c; while ((c = getchar()) != '\n' && c != EOF); }
static int read_int_prompt(const char *prompt) { int v; printf("%s", prompt); while (scanf("%d", &v) != 1) { read_line_discard(); printf("Invalid. %s", prompt); } read_line_discard(); return v; }
static double read_double_prompt(const char *prompt) { double v; printf("%s", prompt); while (scanf("%lf", &v) != 1) { read_line_discard(); printf("Invalid. %s", prompt); } read_line_discard(); return v; }
static void read_string_prompt(const char *prompt, char *buf, size_t n) { printf("%s", prompt); if (!fgets(buf, (int)n, stdin)) { buf[0] = '\0'; return; } size_t L = strlen(buf); if (L > 0 && buf[L-1] == '\n') buf[L-1] = '\0'; }

int main(void) {
    Image *img = NULL; Image *work = NULL; char filename[512];
    while (1) {
        print_menu();
        int opt; if (scanf("%d", &opt) != 1) { read_line_discard(); continue; } read_line_discard();
        if (opt == 0) break;
        if (opt == 1) {
            read_string_prompt("Enter file path to load (PGM P2/P5): ", filename, sizeof(filename));
            Image *tmp = load_pgm(filename);
            if (!tmp) printf("Failed to load %s\n", filename);
            else {
                if (img) free_image(img); if (work) free_image(work);
                img = tmp; work = duplicate_image(img);
                printf("Loaded %s (%dx%d) maxval=%d\n", filename, img->width, img->height, img->maxval);
            }
        } else if (opt == 2) {
            if (!work) { printf("No image loaded.\n"); continue; }
            int factor = read_int_prompt("Enter integer factor (e.g., 2, 3, 4...): ");
            if (factor < 1) { printf("Factor must be >=1\n"); continue; }
            printf("Zoom (z) or Shrink (s)? (z/s): ");
            int c = getchar(); read_line_discard();
            if (c == 'z' || c == 'Z') {
                Image *out = zoom_bilinear(work, factor);
                free_image(work); work = out;
                printf("Zoomed by %dx -> new size: %dx%d\n", factor, work->width, work->height);
            } else if (c == 's' || c == 'S') {
                if (factor < 2) { printf("Shrink factor must be >=2\n"); continue; }
                Image *out = shrink_avg(work, factor);
                free_image(work); work = out;
                printf("Shrunk by %dx -> new size: %dx%d\n", factor, work->width, work->height);
            } else { printf("Invalid choice.\n"); }
        } else if (opt == 3) {
            if (!work) { printf("No image loaded.\n"); continue; }
            printf("Choose filter: 1=box, 2=gaussian, 3=median : ");
            int f; if (scanf("%d", &f) != 1) { read_line_discard(); continue; } read_line_discard();
            if (f == 1) { int k = read_int_prompt("Kernel size (odd >=1, e.g. 3): "); Image *out = box_filter(work, k); free_image(work); work = out; printf("Applied box filter k=%d\n", k); }
            else if (f == 2) { int k = read_int_prompt("Kernel size (odd >=3, e.g. 5): "); double sigma = read_double_prompt("Sigma (e.g. 1.0): "); Image *out = gaussian_filter(work, k, sigma); free_image(work); work = out; printf("Applied gaussian k=%d sigma=%.2f\n", k, sigma); }
            else if (f == 3) { int k = read_int_prompt("Kernel size (odd >=3, e.g. 3): "); Image *out = median_filter(work, k); free_image(work); work = out; printf("Applied median filter k=%d\n", k); }
            else printf("Invalid filter choice.\n");
        } else if (opt == 4) {
            if (!work) { printf("No image loaded.\n"); continue; }
            printf("Choose edge detector: 1=Canny, 2=Sobel, 3=Prewitt : ");
            int e; if (scanf("%d", &e) != 1) { read_line_discard(); continue; } read_line_discard();
            if (e == 1) { double low = read_double_prompt("Low threshold ratio (e.g. 0.05): "); double high = read_double_prompt("High threshold ratio (e.g. 0.15): "); Image *out = canny_edge(work, low, high); free_image(work); work = out; printf("Applied Canny (low=%.3f high=%.3f)\n", low, high); }
            else if (e == 2) { Image *out = sobel_edge(work); free_image(work); work = out; printf("Applied Sobel\n"); }
            else if (e == 3) { Image *out = prewitt_edge(work); free_image(work); work = out; printf("Applied Prewitt\n"); }
            else printf("Invalid choice.\n");
        } else if (opt == 5) {
            if (!work) { printf("No image loaded.\n"); continue; }
            Image *out = compute_lbp(work); free_image(work); work = out; printf("Computed LBP\n");
        } else if (opt == 6) {
            if (!work) { printf("No image loaded.\n"); continue; }
            char outfn[512]; read_string_prompt("Enter filename to save (e.g. out.pgm): ", outfn, sizeof(outfn));
            printf("Save as ASCII P2 or Binary P5? (p2/p5): ");
            char buf[8]; if (!fgets(buf, sizeof(buf), stdin)) buf[0] = '\0'; for (char *s=buf;*s;s++) if (*s=='\n') *s='\0';
            int ascii = (strcmp(buf, "p2")==0 || strcmp(buf, "P2")==0) ? 1 : 0;
            int ok = save_pgm(work, outfn, ascii);
            if (ok) printf("Saved %s (%s)\n", outfn, ascii ? "P2" : "P5"); else printf("Failed to save %s\n", outfn);
        } else if (opt == 7) {
            if (!work) { printf("No image loaded.\n"); continue; }
            printf("Image info: %dx%d maxval=%d\n", work->width, work->height, work->maxval);
        } else if (opt == 8) {
            if (!work) { printf("No image loaded.\n"); continue; }
            unsigned long long csum = image_checksum(work);
            printf("Dimensions: %dx%d  checksum(sum of pixels) = %llu\n", work->width, work->height, csum);
            print_sample(work);
        } else {
            printf("Unknown option.\n");
        }
    }
    if (img) free_image(img); if (work) free_image(work);
    printf("Exiting.\n");
    return 0;
}
