#include <stdio.h>
#include <opencv2/opencv.hpp>
using namespace cv;

int main(int argc, char** argv )
{
    assert( argc == 2 && "usage: displayImg <Image_Path>\n");

    Mat image = imread( argv[1], 1 );

    assert( image.data && "No image data")

    namedWindow("Original Image", WINDOW_AUTOSIZE );
    imshow("Original Image", image);


    //TODO: In image ho la matrice dei pixel dell'immagine
    

    waitKey(0);
    return 0;
}