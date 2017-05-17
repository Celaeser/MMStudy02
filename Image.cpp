#include "Image.h"
#include "ImageProcess.h"
#include <fstream>
#include <stdexcept>

#pragma pack(1)
struct BITMAP_FILE_HEADER {
    std::uint16_t bfType;
    std::uint32_t bfSize;
    std::uint16_t bfReserved1;
    std::uint16_t bfReserved2;
    std::uint32_t bfOffBits;
};

#pragma pack(1)
struct BITMAP_INFO_HEADER {
    std::uint32_t biSize;
    std::int32_t biWidth;
    std::int32_t biHeight;
    std::uint16_t biPlanes;
    std::uint16_t biBitCount;
    std::uint32_t biCompression;
    std::uint32_t biSizeImageMat;
    std::int32_t biXPelsPerMeter;
    std::int32_t biYPelsPerMeter;
    std::uint32_t biClrUsed;
    std::uint32_t biClrImportant;
};
#pragma pack()



ImageMat::ImageMat() :
        channels(0), width(0), height(0), rawData(nullptr), type(BGR)
		{}

ImageMat ImageMat::createFromBMP(const std::string &inputFileURI) {
    ImageMat imageMat;
    BITMAP_FILE_HEADER header;
    BITMAP_INFO_HEADER info;
    std::ifstream file(inputFileURI, std::ifstream::binary);

	if (file.fail()) {
		return imageMat;
	}

    file.read((char *) &header, sizeof(BITMAP_FILE_HEADER))
            .read((char *) &info, sizeof(BITMAP_INFO_HEADER));

    // Get ImageMat info
	imageMat.width = info.biWidth;
	imageMat.height = info.biHeight > 0 ? info.biHeight : -info.biHeight;
	imageMat.channels = info.biBitCount / 8;
	if (imageMat.channels == 1) {
		imageMat.type = ImageMat::Gray;
	}

    // Get ImageMat data
    int lineSkipCount =
            ((imageMat.width * imageMat.channels + 3) & ~0x03)
                                                 - imageMat.width * imageMat.channels; // Skip filled data

	imageMat.rawData = new ImageMat::Byte[imageMat.width * imageMat.height * imageMat.channels];

    file.seekg(header.bfOffBits, file.beg);

    for (int i = 0; i < imageMat.height; ++i) {
        char *rawDataPos = nullptr;

        if (info.biHeight > 0) {
            rawDataPos = (char *) (imageMat.rawData +
                                   (imageMat.height - 1 - i) * imageMat.width * imageMat.channels);
        } else {
            rawDataPos = (char *) (imageMat.rawData +
                                   i * imageMat.width * imageMat.channels);
        }

        file.read(rawDataPos, imageMat.width * imageMat.channels);
        file.seekg(lineSkipCount, file.cur);
    }

    file.close();
    return imageMat;
}

void ImageMat::doCopy(const ImageMat &ImageMat) {
    if (ImageMat.rawData != nullptr) {
        rawData = new Byte[ImageMat.width * ImageMat.height * ImageMat.channels];
        std::memcpy(rawData, ImageMat.rawData,
                    ImageMat.width * ImageMat.height * ImageMat.channels);
    }
}

ImageMat::ImageMat(const ImageMat &imageMat) :
        channels(imageMat.channels), width(imageMat.width),
        height(imageMat.height), rawData(nullptr), type(imageMat.type) {
    doCopy(imageMat);
}

ImageMat &ImageMat::operator=(const ImageMat &imageMat) {
    delete[] rawData;

    channels = imageMat.channels;
    width = imageMat.width;
    height = imageMat.height;
	type = imageMat.type;
    rawData = nullptr;

    doCopy(imageMat);

    return *this;
}

ImageMat::ImageMat(ImageMat &&imageMat) :
        channels(imageMat.channels), width(imageMat.width),
        height(imageMat.height), rawData(imageMat.rawData), type(imageMat.type) {
	imageMat.rawData = nullptr;
}

ImageMat &ImageMat::operator=(ImageMat &&imageMat) {
    delete[] rawData;

    channels = imageMat.channels;
    width = imageMat.width;
    height = imageMat.height;
    rawData = imageMat.rawData;
	type = imageMat.type;

	imageMat.rawData = nullptr;

    return *this;
}

ImageMat::~ImageMat() {
    delete[] rawData;
}

void ImageMat::saveToBMP(const std::string &outputFileURI) {
    BITMAP_FILE_HEADER header;
    BITMAP_INFO_HEADER info;
    std::ofstream file(outputFileURI, std::ofstream::binary);

	// Save no-BGR image
	if (channels == 3 && type != BGR) {
		ImageMat image = cvtColor(*this, BGR);
		image.saveToBMP(outputFileURI);
		return;
	}

    // Fill header
    int lineSize = (width * channels + 3) & ~0x03;
    header.bfType = 0x4d42;
    header.bfSize = channels == 1 ?
                    1078 + lineSize * height :
                    54 + lineSize * height;
    header.bfReserved1 = 0;
    header.bfReserved2 = 0;
    header.bfOffBits = channels == 1 ? 1078 : 54;
    file.write((const char *) &header, 14);

    // Fill info
    info.biSize = 40;
    info.biWidth = width;
    info.biHeight = height;
    info.biPlanes = 1;
    info.biBitCount = 8 * channels;
    info.biCompression = info.biSizeImageMat = 0;
    info.biXPelsPerMeter = info.biYPelsPerMeter = 0;
    info.biClrUsed = info.biClrImportant = channels == 1 ? 256 : 0;
    file.write((const char *) &info, 40);

    // Fill data
    if (channels == 1) {
        char rgbColor[4] = {0};
        for (int i = 0; i < 256; ++i) {
            rgbColor[0] = rgbColor[1] = rgbColor[2] = i;
            file.write(rgbColor, 4);
        }
    }

    const char *rawDataPos = nullptr;
    char emptyZone[4] = {0};
    for (int i = height - 1; i >= 0; --i) {
        rawDataPos = (const char *) rawData + i * width * channels;
        file.write(rawDataPos, width * channels);
        file.write(emptyZone, lineSize - width * channels);
    }
    file.close();
}

ImageMat::ImageMat(uint32_t width, uint32_t height, uint32_t channels, Type type) : 
	width(width), height(height), channels(channels), type(type) {
    rawData = new Byte[this->width * this->height * this->channels];
}