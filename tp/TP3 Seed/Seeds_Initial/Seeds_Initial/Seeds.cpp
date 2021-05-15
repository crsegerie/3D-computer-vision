// Imagine++ project
// Project:  Seeds
// Author:   Pascal Monasse

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <queue>
#include <iostream>
using namespace Imagine;
using namespace std;

/// Min and max disparities
static const float dMin=-30, dMax=-7;

/// Min NCC for a seed
static const float nccSeed=0.95;

/// Radius of patch for correlation
static const int win=(9-1)/2;
/// To avoid division by 0 for constant patch
static const float EPS=0.1f;

/// A seed
struct Seed {
    Seed(int x0, int y0, int d0, float ncc0)
    : x(x0), y(y0), d(d0), ncc(ncc0) {}
    int x,y, d;
    float ncc;
};

/// Order by NCC
bool operator<(const Seed& s1, const Seed& s2) {
    return (s1.ncc<s2.ncc);
}

/// 4-neighbors
static const int dx[]={+1,  0, -1,  0};
static const int dy[]={ 0, -1,  0, +1};

/// Display disparity map
static void displayDisp(const Image<int> disp, Window W, int subW) {
    Image<Color> im(disp.width(), disp.height());
    for(int j=0; j<disp.height(); j++)
        for(int i=0; i<disp.width(); i++) {
            if(disp(i,j)<dMin || disp(i,j)>dMax)
                im(i,j)=Color(0,255,255);
            else {
                int g = 255*(disp(i,j)-dMin)/(dMax-dMin);
                im(i,j)= Color(g,g,g);
            }
        }
    setActiveWindow(W,subW);
    display(im);
    showWindow(W,subW);
}

/// Show 3D window
static void show3D(const Image<Color> im, const Image<int> disp) {
#ifdef IMAGINE_OPENGL // Imagine++ must have been built with OpenGL support...
    // Intrinsic parameters given by Middlebury website
    const float f=3740;
    const float d0=-200; // Doll images cropped by this amount
    const float zoom=2; // Half-size images, should double measured disparity
    const float B=0.160; // Baseline in m
    FMatrix<float,3,3> K(0.0f);
    K(0,0)=-f/zoom; K(0,2)=disp.width()/2;
    K(1,1)= f/zoom; K(1,2)=disp.height()/2;
    K(2,2)=1.0f;
    K = inverse(K);
    K /= K(2,2);
    std::vector<FloatPoint3> pts;
    std::vector<Color> col;
    for(int j=0; j<disp.height(); j++)
        for(int i=0; i<disp.width(); i++)
            if(dMin<=disp(i,j) && disp(i,j)<=dMax) {
                float z = B*f/(zoom*disp(i,j)+d0);
                FloatPoint3 pt((float)i,(float)j,1.0f);
                pts.push_back(K*pt*z);
                col.push_back(im(i,j));
            }
    Mesh mesh(&pts[0], pts.size(), 0,0,0,0,VERTEX_COLOR);
    mesh.setColors(VERTEX, &col[0]);
    Window W = openWindow3D(512,512,"3D");
    setActiveWindow(W);
    showMesh(mesh);
#else
    std::cout << "No 3D: Imagine++ not built with OpenGL support" << std::endl;
#endif
}

/// Correlation between patches centered on (i1,j1) and (i2,j2). The values
/// m1 or m2 are subtracted from each pixel value.
// m1 and m2 are the mean of the patches in the two images
static float correl(const Image<byte>& im1, int i1,int j1,float m1,
                    const Image<byte>& im2, int i2,int j2,float m2)
{
    float dist=0.0f;

    // ------------- TODO -------------
    // the patch must be fully contained in the two images
    if ((i1-win<0)or(j1-win<0)or(i1+win>=im1.width())or(j1+win>=im1.height()))
        cerr<<"patch not fully contained in image 1"<<endl;
    if ((i2-win<0)or(j2-win<0)or(i2+win>=im2.width())or(j2+win>=im2.height()))
        cerr<<"patch not fully contained in image 2"<<endl;


    float var1=0.0f;
    float var2=0.0f;

    // The patch is a square of side 2win+1
    for(int i=-win; i<win+1; i++){
        for(int j=-win; j<win+1; j++){
           dist+=(im1(i+i1,j+j1)-m1)*(im2(i+i2,j+j2)-m2);
           var1+=pow(im1(i+i1,j+j1)-m1, 2);
           var2+=pow(im2(i+i2,j+j2)-m2, 2);
        }
    }
    return dist/(EPS+sqrt(var1*var2));
}

/// Sum of pixel values in patch centered on (i,j).
static float sum(const Image<byte>& im, int i, int j)
{
    if ((i-win<0)or(j-win<0)or(i+win>=im.width())or(j+win>=im.height()))
        cerr<<"patch not fully contained in image"<<endl;

    float s=0.0f;
    // ------------- TODO -------------
    for(int i_=i-win; i_<i+win+1; i_++){
        for(int j_=j-win; j_<j+win+1; j_++){
            s+=im(i_,j_);
        }
    }
    return s;
}

/// Centered correlation of patches of size 2*win+1.
static float ccorrel(const Image<byte>& im1,int i1,int j1,
                     const Image<byte>& im2,int i2,int j2)
{
    float m1 = sum(im1,i1,j1);
    float m2 = sum(im2,i2,j2);
    int w = 2*win+1;
    return correl(im1,i1,j1,m1/(w*w), im2,i2,j2,m2/(w*w));
}

/// Compute disparity map from im1 to im2, but only at points where NCC is
/// above nccSeed. Set to true the seeds and put them in Q.
static void find_seeds(Image<byte> im1, Image<byte> im2,
                       float nccSeed,
                       Image<int>& disp, Image<bool>& seeds,
                       std::priority_queue<Seed>& Q) {
    disp.fill(dMin-1); // so we infer dMin and dMax are authorized disparities
    seeds.fill(false);
    while(! Q.empty())
        Q.pop();

    const int maxy = std::min(im1.height(),im2.height());
    const int refreshStep = (maxy-2*win)*5/100;
    for(int y=win; y+win<maxy; y++) {
        if((y-win-1)/refreshStep != (y-win)/refreshStep)
            std::cout << "Seeds: " << 5*(y-win)/refreshStep <<"%\r"<<std::flush;
        for(int x=win; x+win<im1.width(); x++) {
            // ------------- TODO -------------
            // Hint: just ignore windows that are not fully in image : thanks ! It is so much easier !
            int best_d = dMin-1;
            float best_cc = -1.0f;
            for(int d=dMin; (d<dMax+1)and(x+d+win<im2.width())and(x+d-win>=0); d++) {
                float cc = ccorrel(im1,x,y,im2,x+d,y);
                if (cc>best_cc){
                    best_cc=cc;
                    best_d=d;
                }
            }
            if (best_cc>nccSeed){
                disp(x,y)=best_d;
                seeds(x,y)=true;
                Seed s(x,y,best_d,best_cc);
                Q.push(s);
            }

        }
    }
    std::cout << std::endl;
}

/// Propagate seeds
static void propagate(Image<byte> im1, Image<byte> im2,
                      Image<int>& disp, Image<bool>& seeds,
                      std::priority_queue<Seed>& Q) {
    while(! Q.empty()) {
        Seed s=Q.top();
        Q.pop();
        for(int i=0; i<4; i++) {
            float x=s.x+dx[i], y=s.y+dy[i];
            if(0<=x-win && 0<=y-win &&
               x+win<im2.width() && y+win<im2.height() && // We want the patches associated with the neighbourg pixel fully contained in the images
               ! seeds(x,y)) {
                // ------------- TODO -------------
                // for each adjacent pixel not already in Q, we take the best ccorrelation
                int best_d = s.d; // prior de regularité : Dans l'ancien code int best_d = dMin-1;
                float best_cc = -1.0f;
                for(int d=s.d-1; (d<s.d+2); d++) {
                    if ((x+d+win<im2.width())and(x+d-win>=0)){
                        float cc = ccorrel(im1,x,y,im2,x+d,y);
                        if (cc>best_cc){
                            best_cc=cc;
                            best_d=d;
                        }
                    }

                }
                if (best_d>dMax)
                    best_d=dMax;
                if (best_d<dMin)
                    best_d=dMin;

                disp(x,y)=best_d; // best_d is at least better (or equal) than ccorrel(im1,x,y,im2,x+s.d,y)
                seeds(x,y)=true;
                Seed s(x,y,best_d,best_cc);
                Q.push(s);

            }
        }
    }
}

int main()
{
    // Load and display images
    Image<Color> I1, I2;
    if( ! load(I1, srcPath("im1.jpg")) ||
        ! load(I2, srcPath("im2.jpg")) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }

    if (win >= min(min(I1.width(),I1.height()),min(I2.width(),I2.height())))
        cerr << "At least one of your images are too small or your window is too large" << endl;
        // This condition is necessary in order to write :
        // "for(int i_=max(i-win,0); i_<min(i+win,im.width()); i_++)"

    std::string names[5]={"image 1","image 2","dense","seeds","propagation"};
    Window W = openComplexWindow(I1.width(), I1.height(), "Seeds propagation",
                                 5, names);
    setActiveWindow(W,0);
    display(I1,0,0);
    setActiveWindow(W,1);
    display(I2,0,0);

    Image<int> disp(I1.width(), I1.height());
    Image<bool> seeds(I1.width(), I1.height());
    std::priority_queue<Seed> Q;

    // Dense disparity
    find_seeds(I1, I2, -1.0f, disp, seeds, Q);
    displayDisp(disp,W,2);

    // Only seeds
    find_seeds(I1, I2, nccSeed, disp, seeds, Q);
    displayDisp(disp,W,3);

    // Propagation of seeds
    propagate(I1, I2, disp, seeds, Q);
    displayDisp(disp,W,4);

    // Show 3D (use shift click to animate)
    show3D(I1,disp);

    endGraphics();
    return 0;
}
