
#ifndef IMAGEPROCESS_H
#define IMAGEPROCESS_H

#include <vector>
#include "Image.h"

template<typename T>
struct SquareMatrix {

    std::vector<T> data;
    uint32_t colsOrRows;

    SquareMatrix() {
    }

    SquareMatrix(std::initializer_list<uint8_t> list, uint32_t n) : data(list), colsOrRows(n) {
    }

};

// Predefined DitherMatrix
extern const SquareMatrix<uint8_t> DM2;
extern const SquareMatrix<uint8_t> DM3;
extern const SquareMatrix<uint8_t> DM4;
extern const SquareMatrix<uint8_t> DM8;

/**
 * Calculate the dither matrix which side length equals 2^index
 *
 * @param index     SideLength = 2^index
 * @return          The pointer to dither matrix, <b> Don't forget to delete! </b>
 */
SquareMatrix<uint8_t> *calcDitherMatrix(uint32_t index);

/**
 * Do dithering
 *
 * @param origin            The Origin Image
 * @param ditheringMatrix   The Dithering Matrix
 * @return Dithered image
 */
ImageMat dithering(const ImageMat &origin, const SquareMatrix<uint8_t> *ditheringMatrix = &DM4);

/**
 * Do ordered dithering
 *
 * @param origin                The Origin Image
 * @param ditheringMatrix       The Dithering Matrix
 * @return OrderedDithered image
 */
ImageMat ordered_dithering(const ImageMat &origin, const SquareMatrix<uint8_t> *ditheringMatrix = &DM4);

/**
 * Do dithering
 *
 * @param origin                The Origin Image
 * @param ditheringMatrix       Dithering Matrix in linear-array format
 * @param matrixColsOrRows      Dithering Matrix's cols or rows (side length)
 * @return
 */
ImageMat dithering(const ImageMat &origin, uint8_t *ditheringMatrix, uint32_t matrixColsOrRows);

/**
 *
 * @param origin             The Origin Image
 * @param ditheringMatrix    Dithering Matrix in linear-array format
 * @param colsOrRows         Dithering Matrix's cols or rows (side length)
 * @return
 */
ImageMat ordered_dithering(const ImageMat &origin, uint8_t *ditheringMatrix, uint32_t colsOrRows);

/**
*
* @param origin             The Origin Image
* @param ditheringMatrix    Dithering Matrix in linear-array format
* @param colsOrRows         Dithering Matrix's cols or rows (side length)
* @return
*/
ImageMat ordered_dithering(const ImageMat &origin, uint8_t *ditheringMatrix, uint32_t colsOrRows);

/**
*
* @param origin             The Origin Image
* @param histogram			The histogram text output stream( include the original data and the data after histogram equalization )
* @return					As a gray scale image, there are one image : 
*								-The image after histogram equalization
*							As a colorful image, there are four images( by order ) :
								-The image after histogram equalization in luminance
*								-The chrominance image 1( U for YUV, Cb for YCbCr, H for HSI )
*								-The chrominance image 2( V for YUV, Cr for YCbCr, S for HSI )
*								-The original luminance image
*								-The luminance image after histogram equalization
*/
std::vector<ImageMat> histogram_equalization(const ImageMat &origin, std::ostream& histogram);

/**
* Convert color (origin and output can be the same)
*
* @param origin
* @param outputType
*/
ImageMat cvtColor(const ImageMat &origin, ImageMat::Type outputType);

#endif //IMAGEPROCESS_H