// Imagine++ project
// Project:  Panorama
// Author:   Pascal Monasse
// Date:     2013/10/08

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <sstream>
using namespace Imagine;
using namespace std;

// Record clicks in two images, until right button click
void getClicks(Window w1, Window w2,
               vector<IntPoint2>& pts1, vector<IntPoint2>& pts2) {
    // ------------- TODO/A completer ----------
    int x,y;
    int button=1;
    cout << "Click on a minimum of 4 points on the first image (the one with the car on the right). Click right on the last point."<<endl;
    setActiveWindow(w1);
    while (button != 3) {
        button = getMouse(x,y);
        cout << "Continue."<<endl;
        IntPoint2 pt(x,y);
        pts1.push_back(pt);
    }

    button=1;
    setActiveWindow(w2);
    cout << "Click on the same points on the second image (the one with the car on the left). Click right on the last point."<<endl;
    while (button != 3) {
        button = getMouse(x,y);
        cout << "Continue."<<endl;
        IntPoint2 pt(x,y);
        pts2.push_back(pt);
    }

}

// Return homography compatible with point matches
Matrix<float> getHomography(const vector<IntPoint2>& pts1,
                            const vector<IntPoint2>& pts2) {
    size_t n = min(pts1.size(), pts2.size());
    if(n<4) {
        cout << "Not enough correspondences: " << n << endl;
        return Matrix<float>::Identity(3);
    }
    Matrix<double> A(2*n,8);
    Vector<double> B(2*n);
    // ------------- TODO/A completer ----------
    for(size_t i=0; i<n; i++) {
        A(2*i,0)=pts1[i].x();
        A(2*i,1)=pts1[i].y();
        A(2*i,2)=1;
        A(2*i,3)=0;
        A(2*i,4)=0;
        A(2*i,5)=0;
        A(2*i,6)=-pts2[i].x()*pts1[i].x();
        A(2*i,7)=-pts2[i].x()*pts1[i].y();

        A((2*i)+1,0)=0;
        A((2*i)+1,1)=0;
        A((2*i)+1,2)=0;
        A((2*i)+1,3)=pts1[i].x();
        A((2*i)+1,4)=pts1[i].y();
        A((2*i)+1,5)=1;
        A((2*i)+1,6)=-pts2[i].y()*pts1[i].x();
        A((2*i)+1,7)=-pts2[i].y()*pts1[i].y();

        B[(2*i)]=pts2[i].x();
        B[(2*i)+1]=pts2[i].y();
    }

    B = linSolve(A, B);
    Matrix<float> H(3, 3);
    H(0,0)=B[0]; H(0,1)=B[1]; H(0,2)=B[2];
    H(1,0)=B[3]; H(1,1)=B[4]; H(1,2)=B[5];
    H(2,0)=B[6]; H(2,1)=B[7]; H(2,2)=1;

    // Sanity check
    for(size_t i=0; i<n; i++) {
        float v1[]={(float)pts1[i].x(), (float)pts1[i].y(), 1.0f};
        float v2[]={(float)pts2[i].x(), (float)pts2[i].y(), 1.0f};
        Vector<float> x1(v1,3);
        Vector<float> x2(v2,3);
        x1 = H*x1;
        cout << "sanity check : should be close enouth to zero : "<<endl;
        cout << x1[1]*x2[2]-x1[2]*x2[1] << ' '
             << x1[2]*x2[0]-x1[0]*x2[2] << ' '
             << x1[0]*x2[1]-x1[1]*x2[0] << endl;
    }
    return H;
}

// Grow rectangle of corners (x0,y0) and (x1,y1) to include (x,y)
void growTo(float& x0, float& y0, float& x1, float& y1, float x, float y) {
    if(x<x0) x0=x;
    if(x>x1) x1=x;
    if(y<y0) y0=y;
    if(y>y1) y1=y;    
}

// Panorama construction
void panorama(const Image<Color,2>& I1, const Image<Color,2>& I2,
              Matrix<float> H) {
    Vector<float> v(3);
    float x0=0, y0=0, x1=I2.width(), y1=I2.height();

    v[0]=0; v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=I1.height(); v[2]=1;
    v=H*v;
    cout << "v[2]="<<v[2] <<endl;
    v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=0; v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    // We obtain the coordinate which contain all the projection of the initial image
    cout << "x0 x1 y0 y1=" << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1<<endl;

    Image<Color> I(int(x1-x0), int(y1-y0));
    setActiveWindow( openWindow(I.width(), I.height()) );
    I.fill(WHITE);
    // ------------- TODO/A completer ----------
    for(int i=0; i<I2.width(); i++) {
        for(int j=0; j<I2.height(); j++) {
        I(i-x0,j-y0).r()=I2(i,j).r()/2;
        I(i-x0,j-y0).g()=I2(i,j).g()/2;
        I(i-x0,j-y0).b()=I2(i,j).b()/2;
    }
    }
    display(I,0,0);
    anyClick();
    // push
    if(0){
        for(int i=0; i<I1.width(); i++) {
            for(int j=0; j<I1.height(); j++) {
            v[0]=i; v[1]=j; v[2]=1;
            v=H*v; v/=v[2];
            int a,b;a=(int) v[0]-x0;b=(int) v[1]-y0;
            I(a,b).r()+=I1(i,j).r()/2;
            I(a,b).g()+=I1(i,j).g()/2;
            I(a,b).b()+=I1(i,j).b()/2;
        }
    }}

    // Pull
    if (1){
        Matrix<float> invH(3,3);
        invH = inverse(H);
        for(int i=0; i<I.width(); i++) {
            for(int j=0; j<I.height(); j++) {
                v[0]=i+x0; v[1]=j+y0; v[2]=1;
                v=invH*v; v/=v[2];
                if (((v[0]>0) and (v[0]<I1.width()))and((v[1]>0) and (v[1]<I1.height()))){
                    I(i,j).r()+=I1(v[0],v[1]).r()/2;
                    I(i,j).g()+=I1(v[0],v[1]).g()/2;
                    I(i,j).b()+=I1(v[0],v[1]).b()/2;
                }
            }
        }
    }


    display(I,0,0);
}

// Main function
int main(int argc, char* argv[]) {
    const char* s1 = argc>1? argv[1]: srcPath("image0006.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("image0007.jpg");

    // Load and display images
    Image<Color> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load the images" << endl;
        return 1;
    }
    Window w1 = openWindow(I1.width(), I1.height(), s1);
    display(I1,0,0);
    Window w2 = openWindow(I2.width(), I2.height(), s2);
    setActiveWindow(w2);
    display(I2,0,0);

    // Get user's clicks in images
    vector<IntPoint2> pts1, pts2;
    getClicks(w1, w2, pts1, pts2);

    // Bonus : position of the clicks
    for(int i=0; i<pts2.size(); i++)
        drawCircle(IntPoint2(pts2[i].x(),pts2[i].y()),5,RED);
    setActiveWindow(w1);
    for(int i=0; i<pts1.size(); i++)
        drawCircle(IntPoint2(pts1[i].x(),pts1[i].y()),5,RED);

    vector<IntPoint2>::const_iterator it;
    cout << "pts1="<<endl;
    for(it=pts1.begin(); it != pts1.end(); it++)
        cout << *it << endl;
    cout << "pts2="<<endl;
    for(it=pts2.begin(); it != pts2.end(); it++)
        cout << *it << endl;

    // Compute homography
    Matrix<float> H = getHomography(pts1, pts2);
    cout << "H=" << H/H(2,2);

    // Apply homography
    panorama(I1, I2, H);
    anyClick();

    endGraphics();
    return 0;
}
