#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include "Image.h"
#include "ImageProcess.h"
#include "getopt.h"

#include <array>

#define USAGE "=== Help Info ===\n"\
			  "-h : Show help info\n"\
			  "-i path : Set image file path\n"\
			  "-o path : Set path of output file after histogram equalization( default for outputImageX.bmp )\n"\
			  "-t path : Set path of output text file of histogram( default for histogramX.txt )\n"\
			  "-c color type(YUV \\ YCbCr \\ YIQ \\ HSI) : Set color type( default for YUV )\n"

int main(int argc, char* argv[])
{
	::opterr = 0;
	char optch = 0;
	bool over = false;

	ImageMat::Type colorType = ImageMat::Gray;
	std::string inputPath;
	std::ostringstream outpoutPath;
	std::ostringstream textPath;
	int imageCount = 0;
	std::ifstream fileTester;

	if (argc < 2) {
		std::cout << USAGE << std::endl;
		return 0;
	}

	do {
		fileTester.close();
		outpoutPath.str(std::string());
		outpoutPath << "outputImage" << imageCount << ".bmp";
		++imageCount;
		fileTester.open(outpoutPath.str());
	} while (!fileTester.fail());
	textPath << "histogram" << imageCount - 1 << ".txt";

	while ((optch = getopt(argc, argv, "hi:o:t:c:")) != -1 && !over) {
		switch (optch) {
		case 'h':
		case '?':
			std::cout << USAGE << std::endl;
			over = true;
			break;

		case 'i':
			inputPath = optarg;
			break;

		case 'o':
			outpoutPath.str(optarg);
			break;

		case 't':
			textPath.str(optarg);
			break;

		case 'c':
			if (std::strcmp(optarg, "YUV") == 0) {
				colorType = ImageMat::YUV;
			}
			else if (std::strcmp(optarg, "YCbCr") == 0) {
				colorType = ImageMat::YCbCr;
			}
			else if (std::strcmp(optarg, "YIQ") == 0) {
				colorType = ImageMat::YIQ;
			}
			else if (std::strcmp(optarg, "HSI") == 0) {
				colorType = ImageMat::HSI;
			}
			else{
				std::cerr << "Can't do histogram equalization for type " << optarg << std::endl;
				over = true;
			}
			break;
		}
	}

	if (!over) {
		std::cout << "Reading image from bmp file..." << std::endl;
		ImageMat image = ImageMat::createFromBMP(inputPath);

		std::cout << "Creating histogram text file..." << std::endl;
		std::ofstream histogramText(textPath.str());

		if (image.getRawData() == nullptr) {
			std::cout << "Can't read bmp file " << inputPath << std::endl;
			return 1;
		}

		if (image.getChannels() == 1) {
			if (colorType != ImageMat::Gray) {
				std::cerr << "Can't set color type for gray scale image." << std::endl;
				return 1;
			}
		}
		else {
			std::cout << "Color type changing..." << std::endl;
			if (colorType == ImageMat::Gray) {
				colorType = ImageMat::YUV;
			}

			image = cvtColor(image, colorType);
		}

		std::cout << "Histogram equalization processing..." << std::endl;
		auto imageSet = histogram_equalization(image, histogramText);

		std::cout << "Saving image file..." << std::endl;
		imageSet[0].saveToBMP(outpoutPath.str());

		if (image.getChannels() != 1) {
			const std::string tempPath = outpoutPath.str();
			std::string path;
			path = tempPath.substr(0, tempPath.length() - 4) + "-c0.bmp";
			imageSet[1].saveToBMP(path);
			path = tempPath.substr(0, tempPath.length() - 4) + "-c1.bmp";
			imageSet[2].saveToBMP(path);
			path = tempPath.substr(0, tempPath.length() - 4) + "-lb.bmp";
			imageSet[3].saveToBMP(path);
			path = tempPath.substr(0, tempPath.length() - 4) + "-la.bmp";
			imageSet[4].saveToBMP(path);
		}
	}

	return 0;
}