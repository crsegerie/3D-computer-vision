// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse
// Date:     2013/10/08

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}
vector<Match> randomBatch (vector<Match>& matches, size_t batchSize){
    if (matches.size() < batchSize){
        cerr <<"The SIFT algorithm did not obtain enouth matches " <<endl;
    }
    vector<int> randomIndices;
    bool allDistinct;
    for (size_t i=0;i<batchSize;i++){
        int n;
        do{
            allDistinct = 1;
            int tryN(rand()%matches.size());

            // verifying that n is distinct from the previous selected indices
            for (size_t j=0;j<randomIndices.size();j++){
                if (tryN==randomIndices[j])
                    allDistinct = 0;
            }
            n=tryN;

        }while(allDistinct==0);
        randomIndices.push_back(n);

    }
    vector<Match> randomMatches;
    for (size_t i=0;i<batchSize;i++){
        randomMatches.push_back(matches[randomIndices[i]]);
    }
    return randomMatches;
}

FMatrix<float,3,3> routineGetF (const vector<Match>& matches){
    if(matches.size()<8){
        cerr << "We did not received enouth matches to calculate F" << endl;
    }

    // Normalisation matrix
    FMatrix<float,3,3> N(0.);
    N(0,0) = 0.001;
    N(1,1) = 0.001;
    N(2,2) = 1;

    // Creation of matrix A. Each line of A represents a contraint of a pair of points
    // Even if there is only 8 points, we add a blanc line in A in order to make easily the SVD decomposition
    FMatrix<float,9,9> A(0.);
    for (int i=0;i<8;i++){
        Match match = matches[i];
        DoublePoint3 point1(match.x1,match.y1,1);
        DoublePoint3 point2(match.x2,match.y2,1);

        // Normalizing
        DoublePoint3 normalisedPoint1(N*point1);
        DoublePoint3 normalisedPoint2(N*point2);

        // Filling one line of A
        for (int one=0;one<3;one++){
            for (int two=0;two<3;two++){
                A(i,(3*one)+two) = normalisedPoint1[one]*normalisedPoint2[two];
            }
        }
    }

    // We now extract the smallest eigen value of At*A
    FVector<float,9> S;
    FMatrix<float,9,9> U, Vt;
    svd(A,U,S,Vt);
    FMatrix<float,9,9> V(transpose(Vt));
    FMatrix<float,3,3> F;
    for (int l=0;l<3;l++){
        for (int c=0;c<3;c++){
            F(l,c) = V((3*l)+c,8);
        }
    }

    // Projecting F into the space of rank 2 matrices
    FVector<float,3> S2;
    FMatrix<float,3,3> U2, Vt2;
    svd(F,U2,S2,Vt2);
    S2[2] = 0;
    F = U2 * Diagonal(S2) * Vt2;

    // returning the unormalized F
    return N*F*N;
}

vector<int> computeNumberInliers(const vector<Match>& matches, const FMatrix<float,3,3>& F, const float& distMax){
    vector<int> trueMatches;
    for(size_t i=0; i<matches.size(); i++){
        Match match = matches[i];
        DoublePoint3 p1(match.x1,match.y1,1);
        DoublePoint3 p2(match.x2,match.y2,1);
        FVector<float,3> Ftx(transpose(F)*p1);
        float num = abs(Ftx[0]*p2.x()+Ftx[1]*p2.y()+Ftx[2]);
        float den = sqrt((Ftx[0]*Ftx[0]) + (Ftx[1]*Ftx[1]));
        if (num/den < distMax)
            trueMatches.push_back(i);
    }
    return trueMatches;
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    const size_t batchSize = 8;
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    int bestm=50; // We do not begin at 0 in order to avoid division by 0
    vector<int> bestInliers;
    // --------------- TODO ------------
    cout << "Computing RANSAC algorithm..."<<endl;
    for(int i=0; i<Niter; i++){
        cout << i <<"/"<<Niter<<endl;
        // We take  k = 8 points randomly within the matches vector
        vector<Match> randomMatches;
        randomMatches = randomBatch (matches, batchSize);

        // Calculating the fondamental matrix associated with those 8 points
        FMatrix<float,3,3> F;
        F = routineGetF(randomMatches);

        // Calculating the number m of inliers according to this matrix
        vector<int> inliers = computeNumberInliers(matches, F, distMax);
        int m = inliers.size();

        // Updating the number of needed iterations and the model in case of new ameliorations
        if (m> bestm){
            bestm = m;
            Niter = int(log(BETA)/(log(1-pow((double) m/matches.size(),batchSize))));
            cout << "New best model with "<< bestm <<" inliers found." << endl;
            cout << "New number of needed iterations is : "<<Niter<< endl;
            bestInliers.clear();
            bestInliers = inliers;
            bestF = F;
        }
    }

    cout << "End of Ransac algorithm."<< endl;

    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);


    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    cout <<"Click in one image to show the corresponding epipolar line in other image. Right click to terminate"<<endl;
    while(true) {

        int x,y;
        if(getMouse(x,y) == 3)
            break;
        // --------------- TODO ------------

        // We distinguish the cases where the clicked point is in the left or right image.
        int w = I1.width();
        if (x > w){
            // The user clicked on the right image
            // We draw the epipolar line on the left image
            DoublePoint3 p_prime(x-w, y, 1);
            DoublePoint3 Fp_prime(F*p_prime);

            // the epipolar lines follows the equation x1*a+y1*b+c = 0 with a b c the coeffiscients of FP_prime
            if (Fp_prime[1])
                drawLine(0,-Fp_prime[2]/Fp_prime[1], w,-(Fp_prime[2]+(w*Fp_prime[0]))/Fp_prime[1], RED);
            else
                cerr << "Not implemented ! Bad luck ! The epipolar line is vertical, but try again with a different random seed and it should be ok."<<endl;
        }
        else{
            // The user clicked on the left image
            // We draw the epipolar line on the right image
            DoublePoint3 p(x, y, 1);
            DoublePoint3 ptF(transpose(F)*p); // (pt*F)t = Ft * p
            int w2 = I2.width();
            // the epipolar lines follows the equation x2*a+y2*b+c = 0 with a b c the coeffiscients of FP_prime
            if (ptF[1])
                drawLine(w,-ptF[2]/ptF[1], w+w2,-(ptF[2]+(w2*ptF[0]))/ptF[1], RED);
            else
                cerr << "Not implemented ! Bad luck ! The epipolar line is vertical, but try again with a different random seed and it should be ok."<<endl;
        }

    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    cout << " matches: " << matches.size() << endl;
    cout<<"Please click to begin the Ransac algorithm..."<<endl;
    click();
    
    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);        
    }

    cout << "Please Click to continue..."<< endl;
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
