
/* decal: surface-gradient decal draping on a three-face cliff.
 *
 * 地形是一张高度场 h(x,y)。一株草贴花最初是平面上的单位正方形 (u,v) in [0,1]^2，
 * 通过 r(u,v) = (x0 + d1*u, y0 + d2*v, h(x0 + d1*u, y0 + d2*v)) 贴合到表面。
 * 贴花的变形（雅可比）矩阵由表面梯度直接推出，见 surface_map()。
 *
 * 三种情形：
 *   1. 地表/缓坡贴花：r 是高度场参数化，J 为 3x2，拉伸 = sqrt(1 + |grad h|^2)。
 *   2. 三面垂直悬崖：高度场在阶跃处梯度无穷大，3x2 映射退化，
 *      改用显式墙面参数化 (u,w)，w 为沿墙向下的弧长，
 *      墙面仍是高度场的分段线性部分，J = [L*t, -H*n3]，推导见 wall_map()。
 *   3. 贴花失活：仅当法线接近水平（n_z < NMIN，即坡度超过 70°）时失活。
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* 小向量/小矩阵工具                                                    */
/* ------------------------------------------------------------------ */

typedef struct { float x, y, z; } V3;

static V3 v3(float x, float y, float z) { V3 r = {x, y, z}; return r; }
static float v3dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 v3cross(V3 a, V3 b) {
    return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
static float v3len(V3 a) { return sqrtf(v3dot(a, a)); }
static V3 v3norm(V3 a) { float l = v3len(a); return l > 1e-12f ? v3(a.x / l, a.y / l, a.z / l) : v3(0, 0, 0); }

/* 3x2 变形（雅可比）矩阵，列 c1、c2 分别是 r 对 u、v 的偏导。 */
typedef struct { V3 c1, c2; } Mat32;

/* 右柯西-格林张量 C = J^T J，贴花度量（2x2 对称）。 */
typedef struct { float a, b, c; } Metric22; /* [[a,b],[b,c]] */

static Metric22 metric_of(const Mat32* J) {
    Metric22 m;
    m.a = v3dot(J->c1, J->c1);
    m.b = v3dot(J->c1, J->c2);
    m.c = v3dot(J->c2, J->c2);
    return m;
}

/* 2x2 对称阵特征值（封闭解），输出 lambda1 >= lambda2 > 0。 */
static void metric_eigen(const Metric22* m, float* lam1, float* lam2) {
    float tr = m->a + m->c, det = m->a * m->c - m->b * m->b;
    float d = sqrtf(fmaxf(0.0f, 0.25f * tr * tr - det));
    *lam1 = 0.5f * tr + d;
    *lam2 = 0.5f * tr - d;
    if (*lam2 < 1e-12f) *lam2 = 1e-12f;
}

/* 三档贴花拉伸指标：沿 u 拉伸、沿 v 拉伸、面积拉伸。
 * J 是单位参数域 (u,v) in [0,1]^2 到世界的雅可比；贴花在平面上的原参考
 * 尺寸为 (r1,r2)（地表单位贴花 r1=r2=1；墙面 r1=r2=L）。
 * 归一化度量 C' = D^-1 (J^T J) D^-1，D=diag(r1,r2)，其特征值平方根即
 * 相对于贴花原状的主拉伸倍数。 */
typedef struct { float s_u, s_v, s_area; } Stretch;

static Stretch stretch_normalized(const Mat32* J, float r1, float r2) {
    Metric22 m = metric_of(J);
    Stretch s;
    s.s_u = sqrtf(m.a) / r1;
    s.s_v = sqrtf(m.c) / r2;
    /* |c1 x c2| / (r1*r2)，3D 中单位参考贴花的面积比 */
    s.s_area = v3len(v3cross(J->c1, J->c2)) / (r1 * r2);
    return s;
}

/* 归一化主拉伸：C' = D^-1 C D^-1 的特征值平方根，smax >= smin。 */
static void principal_stretch(const Mat32* J, float r1, float r2,
                              float* smax, float* smin) {
    Metric22 m = metric_of(J);
    Metric22 cn = { m.a / (r1 * r1), m.b / (r1 * r2), m.c / (r2 * r2) };
    float l1, l2;
    metric_eigen(&cn, &l1, &l2);
    *smax = sqrtf(l1);
    *smin = sqrtf(l2);
}

/* ------------------------------------------------------------------ */
/* 高度场                                                              */
/* ------------------------------------------------------------------ */

#define IX0 6.0f
#define IX1 10.0f
#define TOP 6.0f       /* 台地顶高（米） */

/* 解析分段高度场。台地区间 x,y in [6,10]，顶高 TOP：
 *   - 北面 y=10、西面 x=6、东面 x=10 是垂直阶跃（真不连续）；
 *   - 南面台地外 y<6 是 45° 斜坡 h=y（连续过渡，无南面悬崖）；
 *   - 其余外部地面 h=0。 */
typedef struct { float cell; } Field;

static float field_eval(const Field* f, float x, float y) {
    (void)f;
    if (x >= IX0 && x <= IX1 && y >= IX0) {
        if (y <= IX1) return TOP;            /* 台地顶面 */
        return 0.0f;                        /* 台地以北平地（y>10，北悬崖外侧） */
    }
    if (y < IX0) return fmaxf(0.0f, y);     /* 南坡 45°，h=y */
    return 0.0f;                            /* 西/东悬崖外侧平地 */
}

/* 解析中心差分梯度（阶跃两侧各自成立；跨阶跃用 field_jump_grad）。 */
static void field_grad(const Field* f, float x, float y, float* gx, float* gy) {
    float e = 0.5f * f->cell;
    *gx = (field_eval(f, x + e, y) - field_eval(f, x - e, y)) / (2.0f * e);
    *gy = (field_eval(f, x, y + e) - field_eval(f, x, y - e)) / (2.0f * e);
}

/* 跨悬崖面的单侧差分：从台地内侧一点向外侧走一小步 e 跨过不连续面。
 * 垂直面的真实梯度是无穷大（高度在零水平距离内落差 6），有限差分给出的
 * 大数值 |g|=H/e 只是下界；失活判定不依赖该数值，而用墙面的显式法线。
 * axis=0 测 x（西/东墙），axis=1 测 y（北墙）。 */
static void field_jump_grad(const Field* f, float x, float y, int axis, int dir,
                            float* gx, float* gy) {
    float e = 0.25f * f->cell;
    float hi = field_eval(f, x, y);
    float lo = field_eval(f, x + (axis == 0 ? dir * e : 0.0f),
                                 y + (axis == 1 ? dir * e : 0.0f));
    *gx = axis == 0 ? (lo - hi) / e : 0.0f;
    *gy = axis == 1 ? (lo - hi) / e : 0.0f;
}

/* 由梯度直接给出表面单位法线 n = (-gx, -gy, 1)/sqrt(1+gx^2+gy^2)。 */
static V3 normal_from_grad(float gx, float gy) {
    return v3norm(v3(-gx, -gy, 1.0f));
}

/* 失活阈值：法线接近水平（n_z 太小）时贴花失活。
 * NMIN=cos(70°)≈0.342，即坡度角 > 70° 的表面上贴花失活。 */
#define NMIN 0.342f

static int decal_active(float nz) { return nz >= NMIN; }

/* ------------------------------------------------------------------ */
/* 1. 高度场贴花：r(u,v) = (x0+d1*u, y0+d2*v, h(...))                  */
/*    由链式法则直接得到 J 的两列：                                     */
/*      r_u = (d1, 0, d1*hx) = d1*(1,0,hx)                             */
/*      r_v = (0, d2, d2*hy) = d2*(0,1,hy)                             */
/* ------------------------------------------------------------------ */

static Mat32 surface_map(float gx, float gy, float d1, float d2) {
    Mat32 J;
    J.c1 = v3(d1, 0.0f, d1 * gx);
    J.c2 = v3(0.0f, d2, d2 * gy);
    return J;
}

/* ------------------------------------------------------------------ */
/* 2. 垂直悬崖墙面贴花：高度场在阶跃处 3x2 映射退化，改为显式墙面参数化。 */
/*    墙面是高度场的分段线性部分，取墙顶切线 t（水平）与墙面法线 n3。     */
/*    贴花在平面上本是 L×L 的正方形（参数域 (u,v) in [0,1]^2）。          */
/*    上墙后沿墙宽仍占 L，沿墙向下要覆盖整面墙高 H：                     */
/*      r(u,v) = rim + L*u*t - H*v*n3                                   */
/*    J = [L*t, -H*n3]。                                                */
/*    拉伸比要以贴花原尺寸 L 为参考：s_u = L/L = 1，s_v = H/L，           */
/*    s_v 就是"沿崖壁往下拉"的倍数。                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    V3 tangent;  /* 沿墙顶方向的单位水平切线 t */
    V3 normal;   /* 墙面外法线 n3（水平） */
    float L;     /* 贴花原边长（米），参考尺寸 */
    float H;     /* 墙高（米），向下方向需要覆盖的世界长度 */
    char face;   /* 'N' 'W' 'E' */
} Wall;

static Mat32 wall_map(const Wall* wl) {
    Mat32 J;
    J.c1 = v3(wl->L * wl->tangent.x, wl->L * wl->tangent.y, wl->L * wl->tangent.z);
    J.c2 = v3(-wl->H * wl->normal.x, -wl->H * wl->normal.y, -wl->H * wl->normal.z);
    return J;
}

/* ------------------------------------------------------------------ */
/* 三档拉伸分级。阈值取主拉伸 smax（= sqrt(lambda_max(C))）。            */
/*   OK    smax <= 1.5                                                 */
/*   WARN  1.5 < smax <= 2.5                                           */
/*   OVER  smax > 2.5（拉伸超过极限，必须标出所在面）                    */
/* ------------------------------------------------------------------ */

#define WARN_AT 1.5f
#define OVER_AT 2.5f

enum { GRADE_OK = 0, GRADE_WARN = 1, GRADE_OVER = 2 };
static const char* GRADE_NAME[3] = { "OK", "WARN", "OVER" };

static int grade_of(float smax) {
    if (smax > OVER_AT) return GRADE_OVER;
    if (smax > WARN_AT) return GRADE_WARN;
    return GRADE_OK;
}

/* ------------------------------------------------------------------ */
/* 测试场景：4x4 的台地（世界坐标 6..10），顶高 TOP=6 米，单元 1 米。     */
/*   北/西/东三面是垂直阶跃悬崖，南面是 45° 斜坡（|grad|=1）。           */
/* ------------------------------------------------------------------ */

#define CLIFF_H TOP   /* 三面墙的墙高等于台地顶高 */

typedef struct {
    char name[8];
    Mat32 J;
    float r1, r2; /* 贴花原参考尺寸 */
    int grade;
    Stretch st;
    float smax, smin;
    int over;
    char face; /* 0 表示非墙面；否则为 'N' 'W' 'E' */
} DecalReport;

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int feq(float a, float b) { return fabsf(a - b) < 1e-4f; }

static void report_init(DecalReport* r, const char* name, const Mat32* J,
                        float r1, float r2, char face) {
    snprintf(r->name, sizeof(r->name), "%s", name);
    r->J = *J;
    r->r1 = r1;
    r->r2 = r2;
    r->st = stretch_normalized(J, r1, r2);
    principal_stretch(J, r1, r2, &r->smax, &r->smin);
    r->grade = grade_of(r->smax);
    r->over = r->grade == GRADE_OVER;
    r->face = face;
}

static void print_mat(const char* label, const Mat32* J) {
    printf("%-8s J = [(%5.2f %5.2f %5.2f) (%5.2f %5.2f %5.2f)]  (columns dr/du, dr/dv)\n",
           label, J->c1.x, J->c1.y, J->c1.z, J->c2.x, J->c2.y, J->c2.z);
}

/* ------------------------------------------------------------------ */
/* 展开图（surface unfolding net）SVG 生成。                            */
/* 展开规则：各面以真实边长摊平到纸面，铰链边等长。                     */
/*   顶面 TOP 居中；北/西/东三面墙沿各自铰链向外摊开；南坡摊在南侧。     */
/* 每面内画贴花矩形：宽=参考边长 L，高=L*s_v（拉伸后长度）。            */
/* 三档颜色（固定 status token，配符号+斜纹，不只靠颜色）：             */
/*   OK #0ca30c / WARN #fab219 / OVER #d03b3b。                        */
/* ------------------------------------------------------------------ */

#define SVG_S 42.0f          /* 每米像素 */

static const char* GRADE_HEX[3] = { "#0ca30c", "#fab219", "#d03b3b" };
static const char* GRADE_GLYPH[3] = { "✓", "!", "✕" };

/* hatch：非 NULL 时在填充色上叠加 45° 斜纹（颜色之外的冗余编码）。
 * 同名 pattern 只在首次使用时定义。 */
static void svg_rect(FILE* fp, float x, float y, float w, float h,
                     const char* fill, const char* hatch, float sw) {
    if (hatch) {
        const char* line = strcmp(hatch, "hatchWarn") == 0 ? "#8a6200" : "#7a1f1f";
        fprintf(fp, "  <defs><pattern id=\"%s\" width=\"8\" height=\"8\" "
                    "patternUnits=\"userSpaceOnUse\" patternTransform=\"rotate(45)\">"
                    "<rect width=\"8\" height=\"8\" fill=\"%s\"/>"
                    "<line x1=\"0\" y1=\"0\" x2=\"0\" y2=\"8\" stroke=\"%s\" "
                    "stroke-width=\"2\"/></pattern></defs>\n", hatch, fill, line);
    }
    fprintf(fp, "  <rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" "
                "fill=\"%s\" stroke=\"#0b0b0b\" stroke-width=\"%.2f\"/>\n",
            x, y, w, h, hatch ? "none" : fill, sw);
    if (hatch)
        fprintf(fp, "  <rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" "
                    "fill=\"url(#%s)\" stroke=\"#0b0b0b\" stroke-width=\"%.2f\"/>\n",
                x, y, w, h, hatch, sw);
}

static void svg_text(FILE* fp, float x, float y, int size, const char* anchor,
                     const char* fill, const char* txt) {
    fprintf(fp, "  <text x=\"%.1f\" y=\"%.1f\" font-family=\"system-ui,sans-serif\" "
                "font-size=\"%d\" text-anchor=\"%s\" fill=\"%s\">%s</text>\n",
            x, y, size, anchor, fill, txt);
}

/* 单个墙面贴花标签：符号 + 面名 + 拉伸倍数。 */
static void svg_decal_label(FILE* fp, float cx, float cy, const DecalReport* r) {
    char buf[96];
    snprintf(buf, sizeof(buf), "%s face %c  s_v=%.2f",
             GRADE_GLYPH[r->grade], r->face, r->st.s_v);
    svg_text(fp, cx, cy, 15, "middle",
             r->grade == GRADE_WARN ? "#0b0b0b" : "#ffffff", buf);
}

static void write_unfolding_svg(const char* path, const DecalReport* walls,
                                const DecalReport* ground, const DecalReport* slope) {
    FILE* fp = fopen(path, "w");
    if (!fp) { perror(path); return; }

    const float CW = 940.0f, CH = 880.0f;
    float top = 4.0f * SVG_S;                 /* 顶面 4m = 168px */
    float wallH = CLIFF_H * SVG_S;            /* 墙高 6m = 252px */
    float slopeLen = CLIFF_H * 1.414214f * SVG_S; /* 南坡 6√2 m ≈ 356px */
    float tx = 300.0f, ty = 300.0f;           /* 顶面左上角 */

    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" "
                "viewBox=\"0 0 %.0f %.0f\" font-family=\"system-ui,sans-serif\">\n",
                CW, CH, CW, CH);
    fprintf(fp, "  <rect width=\"%.0f\" height=\"%.0f\" fill=\"#fcfcfb\"/>\n", CW, CH);

    /* 标题 */
    svg_text(fp, CW / 2, 34, 22, "middle", "#0b0b0b",
             "悬崖贴花表面展开图（unfolding net，真实边长摊平）");
    svg_text(fp, CW / 2, 58, 13, "middle", "#52514e",
             "铰链边等长；贴花矩形高 = 参考边长 L × 垂直拉伸 s_v；颜色/符号/斜纹三重编码拉伸档位");

    /* ---- 各面摊平区域（先画面，再画贴花）---- */
    float nx = tx, ny = ty - wallH;           /* 北墙：铰链=顶面北缘，向上摊开 */
    svg_rect(fp, nx, ny, top, wallH, "#f0efec", NULL, 1.0f);
    float wx = tx - wallH, wy = ty;           /* 西墙：铰链=西缘，向左摊开 */
    svg_rect(fp, wx, wy, wallH, top, "#f0efec", NULL, 1.0f);
    float ex = tx + top, ey = ty;             /* 东墙：铰链=东缘，向右摊开 */
    svg_rect(fp, ex, ey, wallH, top, "#f0efec", NULL, 1.0f);
    float sx = tx, sy = ty + top;             /* 南坡：铰链=南缘，向下摊开 */
    svg_rect(fp, sx, sy, top, slopeLen, "#f7f6f1", NULL, 1.0f);

    /* ---- 顶面（含一株地面单位贴花，无拉伸）---- */
    svg_rect(fp, tx, ty, top, top, "#e9e7df", NULL, 1.5f);
    float g = SVG_S;
    svg_rect(fp, tx + 0.4f * SVG_S, ty + 0.4f * SVG_S, g, g,
             GRADE_HEX[ground->grade], NULL, 1.0f);
    svg_text(fp, tx + 0.4f * SVG_S + g / 2, ty + 0.4f * SVG_S + g / 2 + 5, 12,
             "middle", "#ffffff", "✓1.0");

    /* ---- 北墙贴花（OK，L=4，s_v=1.5 -> 高 6m=满墙高）---- */
    {
        const DecalReport* r = &walls[0];
        float Lpx = 4.0f * SVG_S, hpx = Lpx * r->st.s_v;
        float dx = nx + (top - Lpx) / 2, dy = ny + wallH - hpx; /* 锚铰链边向下展开 */
        svg_rect(fp, dx, dy, Lpx, hpx, GRADE_HEX[r->grade], NULL, 1.0f);
        svg_decal_label(fp, dx + Lpx / 2, dy + hpx / 2, r);
    }
    /* ---- 西墙贴花（WARN，L=3，s_v=2.0 -> 高 6m），斜纹兜底 ---- */
    {
        const DecalReport* r = &walls[1];
        float Lpx = 3.0f * SVG_S, hpx = Lpx * r->st.s_v;
        float dx = wx + wallH - hpx, dy = wy + (top - Lpx) / 2;
        svg_rect(fp, dx, dy, hpx, Lpx, GRADE_HEX[r->grade], "hatchWarn", 1.0f);
        svg_decal_label(fp, dx + hpx / 2, dy + Lpx / 2, r);
    }
    /* ---- 东墙贴花（OVER，L=2，s_v=3.0 -> 高 6m），斜纹+粗框醒目标注 ---- */
    {
        const DecalReport* r = &walls[2];
        float Lpx = 2.0f * SVG_S, hpx = Lpx * r->st.s_v;
        float dx = ex, dy = ey + (top - Lpx) / 2;
        svg_rect(fp, dx, dy, hpx, Lpx, GRADE_HEX[r->grade], "hatchOver", 3.0f);
        svg_decal_label(fp, dx + hpx / 2, dy + Lpx / 2, r);
    }
    /* ---- 南坡贴花（45°，s_v=√2）---- */
    {
        float Lpx = 2.0f * SVG_S, hpx = Lpx * slope->st.s_v;
        float dx = sx + (top - Lpx) / 2, dy = sy;
        svg_rect(fp, dx, dy, Lpx, hpx, GRADE_HEX[slope->grade], NULL, 1.0f);
        char buf[80];
        snprintf(buf, sizeof(buf), "✓ slope s_v=%.2f", slope->st.s_v);
        svg_text(fp, dx + Lpx / 2, dy + hpx / 2 + 5, 13, "middle", "#ffffff", buf);
    }

    /* ---- 面名与铰链标注 ---- */
    svg_text(fp, tx + top / 2, ty + top / 2 - 8, 15, "middle", "#52514e", "台地顶面 TOP");
    svg_text(fp, tx + top / 2, ty + top / 2 + 12, 12, "middle", "#898781", "4m × 4m");
    svg_text(fp, nx + top / 2, ny + 20, 15, "middle", "#0b0b0b", "北墙 N（铰链 = 北缘）");
    svg_text(fp, wx + 16, wy + top / 2, 13, "middle", "#0b0b0b", "西墙 W");
    svg_text(fp, ex + wallH - 16, ey + top / 2, 13, "middle", "#0b0b0b", "东墙 E");
    svg_text(fp, sx + top / 2, sy + slopeLen - 12, 13, "middle", "#52514e",
             "南坡 45°（坡长 6√2 m，连续面，贴花不失效）");

    /* ---- 超限面警示框（右下角，不与南坡重叠：南坡底 ≈ 824，框放其下）---- */
    float bx = 300.0f, by = 836.0f;
    fprintf(fp, "  <rect x=\"%.1f\" y=\"%.1f\" width=\"360\" height=\"36\" rx=\"8\" "
                "fill=\"#fff5f5\" stroke=\"#d03b3b\" stroke-width=\"2\"/>\n", bx, by);
    svg_text(fp, bx + 14, by + 23, 14, "start", "#d03b3b",
             "✕ OVER 极限（s_v&gt;2.5）→ 所在面：东墙 E，L=2m 盖 6m 墙高，H/L=3.00");

    /* ---- 图例（右侧，x 从 760 起，避开东墙右缘 720）---- */
    float lx = 760.0f, ly = 300.0f;
    svg_text(fp, lx, ly - 20, 15, "start", "#0b0b0b", "拉伸档位（主拉伸 smax）");
    const char* desc[3] = { "OK    smax &lt;= 1.5", "WARN  1.5 &lt; smax &lt;= 2.5",
                            "OVER  smax &gt; 2.5" };
    for (int i = 0; i < 3; i++) {
        float yy = ly + 8 + i * 36;
        fprintf(fp, "  <rect x=\"%.1f\" y=\"%.1f\" width=\"26\" height=\"24\" fill=\"%s\" "
                    "stroke=\"#0b0b0b\" stroke-width=\"%g\"/>\n",
                lx, yy, GRADE_HEX[i], i == 2 ? 2.5 : 1.0);
        if (i >= 1)
            fprintf(fp, "  <line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" "
                        "stroke=\"%s\" stroke-width=\"2\"/>\n",
                        lx + 2, yy + 22, lx + 24, yy + 2, i == 1 ? "#8a6200" : "#7a1f1f");
        svg_text(fp, lx + 13, yy + 18, 13, "middle",
                 i == 1 ? "#0b0b0b" : "#ffffff", GRADE_GLYPH[i]);
        svg_text(fp, lx + 36, yy + 17, 13, "start", "#0b0b0b", desc[i]);
    }
    svg_text(fp, lx, ly + 128, 12, "start", "#898781", "符号 ✓ / ! / ✕ 与斜纹为颜色之外的冗余编码");
    svg_text(fp, lx, ly + 148, 12, "start", "#898781", "失活：n_z &lt; cos70° ≈ 0.342 时贴花不投影");

    fprintf(fp, "</svg>\n");
    fclose(fp);
}

int main(void) {
    Field f = { 1.0f };

    /* ---- 梯度与法线：南坡 |grad|=1，n=(0,-1,1)/sqrt2，n_z≈0.707 ---- */
    float gx, gy;
    field_grad(&f, 3.0f, 3.0f, &gx, &gy);
    check("ramp_grad", feq(gx, 0.0f) && feq(gy, 1.0f));
    V3 nr = normal_from_grad(gx, gy);
    check("ramp_normal", feq(nr.z, 0.707107f) && feq(nr.y, -0.707107f));
    check("ramp_active", decal_active(nr.z) == 1);

    /* 台地顶面梯度为 0，法线竖直，贴花恒激活。 */
    field_grad(&f, 8.0f, 8.0f, &gx, &gy);
    check("flat_grad", feq(gx, 0.0f) && feq(gy, 0.0f));
    check("flat_active", decal_active(normal_from_grad(gx, gy).z) == 1);

    /* 三面垂直悬崖：高度场在墙上不连续，grad h -> 无穷，法线严格水平 n_z=0，
     * 小于 NMIN，贴花在三个面上都失活。有限差分梯度随步长增大（这里 H/e=24），
     * 仅作阶跃存在性的佐证；失活判定用墙面的显式法线。 */
    struct { char face; float x, y; int axis, dir; V3 wn; } cliffs[3] = {
        { 'N', 8.0f, IX1 - 0.1f, 1,  1, v3(0, 1, 0) },  /* 北墙 y=10，朝外 +y 平地 */
        { 'W', IX0 + 0.1f, 8.0f, 0, -1, v3(-1, 0, 0) }, /* 西墙 x=6，朝外 -x */
        { 'E', IX1 - 0.1f, 8.0f, 0,  1, v3(1, 0, 0) },  /* 东墙 x=10，朝外 +x */
    };
    for (int i = 0; i < 3; i++) {
        field_jump_grad(&f, cliffs[i].x, cliffs[i].y, cliffs[i].axis, cliffs[i].dir, &gx, &gy);
        float gmag = sqrtf(gx * gx + gy * gy);
        char nm[24];
        snprintf(nm, sizeof(nm), "cliff%c_jump", cliffs[i].face);
        check(nm, gmag > 1.0f / NMIN); /* 有限差分已超出失活阈值对应的梯度 */
        snprintf(nm, sizeof(nm), "cliff%c_inactive", cliffs[i].face);
        check(nm, cliffs[i].wn.z == 0.0f && decal_active(cliffs[i].wn.z) == 0);
    }

    /* 临界角检查：NMIN 对应 70°，69° 激活、71° 失活。 */
    check("threshold_69_active", decal_active(cosf(69.0f * 3.14159265f / 180.0f)) == 1);
    check("threshold_71_inactive", decal_active(cosf(71.0f * 3.14159265f / 180.0f)) == 0);

    /* ---- 高度场贴花变形矩阵：缓坡 45°，d1=d2=1 ---- */
    Mat32 Jr = surface_map(0.0f, 1.0f, 1.0f, 1.0f);
    check("ramp_J_c1", feq(Jr.c1.x, 1.0f) && feq(Jr.c1.z, 0.0f));
    check("ramp_J_c2", feq(Jr.c2.y, 1.0f) && feq(Jr.c2.z, 1.0f));
    Stretch sr = stretch_normalized(&Jr, 1.0f, 1.0f);
    /* 等高线方向不拉伸；下坡方向 sqrt(1+1)=sqrt2；面积 = sqrt2。 */
    check("ramp_stretch", feq(sr.s_u, 1.0f) && feq(sr.s_v, 1.414214f) && feq(sr.s_area, 1.414214f));
    float rmax, rmin;
    principal_stretch(&Jr, 1.0f, 1.0f, &rmax, &rmin);
    check("ramp_principal", feq(rmax, 1.414214f) && feq(rmin, 1.0f));
    check("ramp_grade_ok", grade_of(rmax) == GRADE_OK);

    /* 平面单位贴花：J 为单位阵，三项指标全 1。 */
    Mat32 Jf = surface_map(0.0f, 0.0f, 1.0f, 1.0f);
    Stretch sf = stretch_normalized(&Jf, 1.0f, 1.0f);
    check("flat_identity", feq(sf.s_u, 1.0f) && feq(sf.s_v, 1.0f) && feq(sf.s_area, 1.0f));

    /* ---- 三面悬崖，各贴一次；墙高 6，贴花外伸宽度 L 不同 => 三档 ---- */
    Wall walls[3];
    memset(walls, 0, sizeof(walls));
    /* 墙顶锚点（注释备查，平移不进雅可比）：
     *   N: rim=(6,10,6)  W: rim=(6,6,6)  E: rim=(10,6,6)。 */
    /* 北墙 y=10：墙顶沿 x 方向，外法线朝 +y。L=4 => s_v=6/4=1.5（OK 档上界）。 */
    walls[0] = (Wall){ v3(1, 0, 0), v3(0, 1, 0), 4.0f, CLIFF_H, 'N' };
    /* 西墙 x=6：墙顶沿 y 方向，外法线朝 -x。L=3 => s_v=6/3=2.0（WARN）。 */
    walls[1] = (Wall){ v3(0, 1, 0), v3(-1, 0, 0), 3.0f, CLIFF_H, 'W' };
    /* 东墙 x=10：墙顶沿 y 方向，外法线朝 +x。L=2 => s_v=6/2=3.0（OVER）。 */
    walls[2] = (Wall){ v3(0, 1, 0), v3(1, 0, 0), 2.0f, CLIFF_H, 'E' };

    DecalReport rep[3];
    memset(rep, 0, sizeof(rep));
    for (int i = 0; i < 3; i++) {
        char nm[8];
        snprintf(nm, sizeof(nm), "wall_%c", walls[i].face);
        Mat32 Jw = wall_map(&walls[i]);
        report_init(&rep[i], nm, &Jw, walls[i].L, walls[i].L, walls[i].face);
    }

    /* 墙面正交：两列点积为 0，沿墙宽方向不拉伸，面积比 = s_u*s_v。 */
    for (int i = 0; i < 3; i++) {
        char nm[20];
        snprintf(nm, sizeof(nm), "wall%c_orthogonal", walls[i].face);
        check(nm, feq(v3dot(rep[i].J.c1, rep[i].J.c2), 0.0f));
    }

    /* 北墙：s_u=1, s_v=1.5, 面积 1.5，OK 档（1.5 不大于阈值）。 */
    check("wallN_stretch", feq(rep[0].st.s_u, 1.0f) && feq(rep[0].st.s_v, 1.5f)
                            && feq(rep[0].st.s_area, 1.5f));
    check("wallN_grade", rep[0].grade == GRADE_OK);

    /* 西墙：s_v=2.0，WARN 档。 */
    check("wallW_stretch", feq(rep[1].st.s_u, 1.0f) && feq(rep[1].st.s_v, 2.0f)
                            && feq(rep[1].st.s_area, 2.0f));
    check("wallW_grade", rep[1].grade == GRADE_WARN);

    /* 东墙：s_v=3.0，OVER 档，必须能报出所在面。 */
    check("wallE_stretch", feq(rep[2].st.s_u, 1.0f) && feq(rep[2].st.s_v, 3.0f)
                            && feq(rep[2].st.s_area, 3.0f));
    check("wallE_grade_over", rep[2].grade == GRADE_OVER && rep[2].over == 1);

    /* ---- 交付物：变形矩阵 ---- */
    printf("\n--- decal deformation matrices (columns = dr/du, dr/dv) ---\n");
    print_mat("ground", &Jf);
    print_mat("slope45", &Jr);
    for (int i = 0; i < 3; i++) {
        char label[16];
        snprintf(label, sizeof(label), "wall_%c", walls[i].face);
        print_mat(label, &rep[i].J);
    }

    /* ---- 交付物：三档拉伸指标表 ---- */
    printf("\n--- stretch metrics: s_u (along u), s_v (along v), s_area ---\n");
    DecalReport all[5];
    report_init(&all[0], "ground", &Jf, 1.0f, 1.0f, 0);
    report_init(&all[1], "slope45", &Jr, 1.0f, 1.0f, 0);
    for (int i = 0; i < 3; i++) all[i + 2] = rep[i];
    printf("%-8s %6s %6s %7s %7s %6s  %s\n",
           "decal", "s_u", "s_v", "s_area", "smax", "smin", "grade");
    for (int i = 0; i < 5; i++) {
        printf("%-8s %6.3f %6.3f %7.3f %7.3f %6.3f  %s%s\n",
               all[i].name, all[i].st.s_u, all[i].st.s_v, all[i].st.s_area,
               all[i].smax, all[i].smin, GRADE_NAME[all[i].grade],
               all[i].over ? "  <-- OVER LIMIT" : "");
        if (all[i].over)
            printf("           \\__ over-limit decal located on face %c\n", all[i].face);
    }

    printf("\nover_limit_faces=");
    int any = 0;
    for (int i = 0; i < 3; i++) if (rep[i].over) { printf("%c ", walls[i].face); any = 1; }
    if (!any) printf("none");
    printf("\n");

    /* ---- 交付物：展开图 SVG ---- */
    const char* svg_path = "unfolding.svg";
    write_unfolding_svg(svg_path, rep, &all[0], &all[1]);
    printf("unfolding_svg=%s\n", svg_path);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
