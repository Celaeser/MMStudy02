#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <array>
#include <utility>
#include <cstdint>
#include "ImageProcess.h"

static std::vector<uint8_t> he_helper(const ImageMat &origin, const int channel, std::ostream& histogram)
{
	std::vector<uint8_t> mapData(256, 0);
	std::array<std::pair<uint32_t, uint32_t>, 256> accuData = { std::make_pair(0, 0) };
	std::array<uint32_t, 256> countData = { 0 };
	const uint32_t TOTAL = origin.getWidth() * origin.getHeight();

	//Count
	for (uint32_t i = 0; i < origin.getWidth(); ++i) {
		for (uint32_t j = 0; j < origin.getHeight(); ++j) {
			++accuData[(size_t)origin.getPixiv(i, j)[channel]].first;
		}
	}

	//Record before equalization
	histogram << "before equalization : \n";
	for (int i = 0; i < 256; ++i) {
		histogram << i << " : " << accuData[i].first << '\n';
	}

	//Accumulate
	accuData[0].second = accuData[0].first;
	for (int i = 1; i < 256; ++i) {
		accuData[i].second = accuData[i - 1].second + accuData[i].first;
	}

	//Equalization
	for (int i = 1; i < 256; ++i) {
		mapData[i] = (uint8_t)std::round(accuData[i].second * 255.0 / TOTAL);
		countData[mapData[i]] = accuData[i].first;
	}

	//Record after equalization
	histogram << "\n\nafter equalization : \n";
	for (int i = 0; i < 256; ++i) {
		histogram << i << " : " << countData[i] << '\n';
	}

	return mapData;
}

std::vector<ImageMat> histogram_equalization(const ImageMat &origin, std::ostream& histogram)
{
	std::vector<ImageMat> imageSet;
	std::vector<uint8_t> mapData;
	const uint32_t TOTAL = origin.getWidth() * origin.getHeight();

	if (origin.getType() == ImageMat::Gray) {
		//Get pixiv map data and output histogram
		mapData = he_helper(origin, 0, histogram);

		//Generate image
		imageSet.emplace_back(origin.getWidth(), origin.getHeight(), 1, ImageMat::Gray);
		for (uint32_t i = 0; i < origin.getWidth(); ++i) {
			for (uint32_t j = 0; j < origin.getHeight(); ++j) {
				imageSet.back().getPixiv(i, j)[0] = mapData[(size_t)origin.getPixiv(i, j)[0]];
			}
		}
	}
	else {
		size_t lumChannel = 0;
		size_t chrChannel[2] = { 0 };

		switch (origin.getType()) {
		case ImageMat::HSI:
			lumChannel = 2;
			chrChannel[0] = 0;
			chrChannel[1] = 1;
			break;

		case ImageMat::YUV:
		case ImageMat::YIQ:
		case ImageMat::YCbCr:
			lumChannel = 0;
			chrChannel[0] = 1;
			chrChannel[1] = 2;
			break;

		default:
			throw std::runtime_error("can't do histogram equalization for images of this type.");
			break;
		}

		//Get pixiv map data and output histogram
		mapData = he_helper(origin, lumChannel, histogram);

		//Generate image
		imageSet.emplace_back(origin.getWidth(), origin.getHeight(), 3, origin.getType());
		for (uint32_t i = 0; i < origin.getWidth(); ++i) {
			for (uint32_t j = 0; j < origin.getHeight(); ++j) {
				imageSet.back().getPixiv(i, j)[lumChannel] = mapData[(size_t)origin.getPixiv(i, j)[lumChannel]];
				imageSet.back().getPixiv(i, j)[chrChannel[0]] = origin.getPixiv(i, j)[chrChannel[0]];
				imageSet.back().getPixiv(i, j)[chrChannel[1]] = origin.getPixiv(i, j)[chrChannel[1]];
			}
		}

		imageSet.emplace_back(origin.getWidth(), origin.getHeight(), 1, ImageMat::Gray);
		for (uint32_t i = 0; i < origin.getWidth(); ++i) {
			for (uint32_t j = 0; j < origin.getHeight(); ++j) {
				imageSet.back().getPixiv(i, j)[0] = origin.getPixiv(i, j)[chrChannel[0]];
			}
		}

		imageSet.emplace_back(origin.getWidth(), origin.getHeight(), 1, ImageMat::Gray);
		for (uint32_t i = 0; i < origin.getWidth(); ++i) {
			for (uint32_t j = 0; j < origin.getHeight(); ++j) {
				imageSet.back().getPixiv(i, j)[0] = origin.getPixiv(i, j)[chrChannel[1]];
			}
		}

		imageSet.emplace_back(origin.getWidth(), origin.getHeight(), 1, ImageMat::Gray);
		for (uint32_t i = 0; i < origin.getWidth(); ++i) {
			for (uint32_t j = 0; j < origin.getHeight(); ++j) {
				imageSet.back().getPixiv(i, j)[0] = origin.getPixiv(i, j)[lumChannel];
			}
		}

		imageSet.emplace_back(origin.getWidth(), origin.getHeight(), 1, ImageMat::Gray);
		for (uint32_t i = 0; i < origin.getWidth(); ++i) {
			for (uint32_t j = 0; j < origin.getHeight(); ++j) {
				imageSet.back().getPixiv(i, j)[0] = mapData[(size_t)origin.getPixiv(i, j)[lumChannel]];
			}
		}
	}

	return imageSet;
}

static bool isSameSize(const ImageMat &origin, const ImageMat &output) {
	if (origin.getWidth() != output.getWidth()) {
		return false;
	}
	if (origin.getHeight() != output.getHeight()) {
		return false;
	}
	return true;
}

static ImageMat cvtColorToBGR(const ImageMat &origin)
{
	ImageMat output = ImageMat(origin.getWidth(), origin.getHeight(), 3, ImageMat::BGR);
	auto checker = [](ImageMat::Byte* ch) {
		for (int i = 0; i < 3; ++i) {
			if (ch[i] > 255.0) ch[i] = 255.0;
			if (ch[i] < 0) ch[i] = 0;
		}
	};

	// Do convert
	uint32_t cols = origin.getWidth();
	uint32_t rows = origin.getHeight();
	const ImageMat::Byte *originData = origin.getRawData();
	ImageMat::Byte *outputData = output.getRawData();

	switch (origin.getType()) {
	case ImageMat::BGR:
		output = origin;
		break;
	case ImageMat::YUV: {
		for (uint32_t i = 0; i < origin.getHeight(); ++i) {
			for (uint32_t j = 0; j < origin.getWidth(); ++j) {
				const ImageMat::Byte *originPix = originData + i * cols * 3 + j * 3;
				ImageMat::Byte *outputPix = outputData + i * cols * 3 + j * 3;
				ImageMat::Byte Y = originPix[0];
				double U = (originPix[1] - 128.0) / 0.5 * 0.437;
				double V = (originPix[2] - 128.0) / 0.5 * 0.615;
				// R
				outputPix[2] = (ImageMat::Byte)(1 * Y + 0.000 * U + 1.140 * V);
				// G
				outputPix[1] = (ImageMat::Byte)(1 * Y - 0.395 * U - 0.581 * V);
				// B
				outputPix[0] = (ImageMat::Byte)(1 * Y + 2.032 * U + 0.000 * V);

				checker(outputPix);
			}
		}
		break;
	}
	case ImageMat::YCbCr: {
		for (uint32_t i = 0; i < origin.getHeight(); ++i) {
			for (uint32_t j = 0; j < origin.getWidth(); ++j) {
				const ImageMat::Byte *originPix = originData + i * cols * 3 + j * 3;
				ImageMat::Byte *outputPix = outputData + i * cols * 3 + j * 3;
				ImageMat::Byte Y = originPix[0];
				ImageMat::Byte Cb = originPix[1];
				ImageMat::Byte Cr = originPix[2];
				// R
				outputPix[2] = (ImageMat::Byte)(1 * Y + 0.00000 * (Cb - 128) + 1.40200 * (Cr - 128));
				// G
				outputPix[1] = (ImageMat::Byte)(1 * Y - 0.34414 * (Cb - 128) - 0.71414 * (Cr - 128));
				// B
				outputPix[0] = (ImageMat::Byte)(1 * Y + 1.77200 * (Cb - 128) + 0.00000 * (Cr - 128));

				checker(outputPix);
			}
		}
		break;
	}
	case ImageMat::YIQ: {
		for (uint32_t i = 0; i < origin.getHeight(); ++i) {
			for (uint32_t j = 0; j < origin.getWidth(); ++j) {
				const ImageMat::Byte *originPix = originData + i * cols * 3 + j * 3;
				ImageMat::Byte *outputPix = outputData + i * cols * 3 + j * 3;
				ImageMat::Byte Y = originPix[0];
				double I = (originPix[1] - 128.0) / 0.5 * 0.596;
				double Q = (originPix[2] - 128.0) / 0.5 * 0.523;
				// R
				outputPix[2] = (ImageMat::Byte)(1 * Y + 0.956 * I + 0.621 * Q);
				// G
				outputPix[1] = (ImageMat::Byte)(1 * Y - 0.272 * I - 0.647 * Q);
				// B
				outputPix[0] = (ImageMat::Byte)(1 * Y - 1.106 * I + 1.703 * Q);

				checker(outputPix);
			}
		}
		break;
	}
	case ImageMat::HSI: {
		for (uint32_t i = 0; i < origin.getHeight(); ++i) {
			for (uint32_t j = 0; j < origin.getWidth(); ++j) {
				// FIXME: some problems here
				const ImageMat::Byte *originPix = originData + i * cols * 3 + j * 3;
				ImageMat::Byte *outputPix = outputData + i * cols * 3 + j * 3;
				double H = originPix[0] / 255.0 * 360;
				double S = originPix[1] / 255.0;
				double I = originPix[2];
				if (H <= 120) {
					// B
					outputPix[0] = (ImageMat::Byte)(I * (1 - S));
					// R
					outputPix[2] = (ImageMat::Byte)(I * (1 + S * cos(M_PI / 180 * H) / cos(M_PI / 180 * (60 - H))));
					// G
					outputPix[1] = (ImageMat::Byte)(3 * I - (outputPix[0] + outputPix[2]));
				}
				else if (H <= 240) {
					H -= 120;
					// R
					outputPix[2] = (ImageMat::Byte)(I * (1 - S));
					// G
					outputPix[1] = (ImageMat::Byte)(I * (1 + S * cos(M_PI / 180 * H) / cos(M_PI / 180 * (60 - H))));
					// B
					outputPix[0] = (ImageMat::Byte)(3 * I - (outputPix[2] + outputPix[1]));

				}
				else {
					H -= 240;
					// G
					outputPix[1] = (ImageMat::Byte)(I * (1 - S));
					// B
					outputPix[0] = (ImageMat::Byte)(I * (1 + S * cos(M_PI / 180 * H) / cos(M_PI / 180 * (60 - H))));
					// R
					outputPix[2] = (ImageMat::Byte)(3 * I - (outputPix[1] + outputPix[0]));
				}

				checker(outputPix);
			}
		}
		break;
	}
	case ImageMat::Gray:
		throw std::runtime_error("Can not convert Gray to BGR");
		break;
	default:
		throw std::runtime_error("Unsupported Type");
		break;
	}
	output.setType(ImageMat::BGR);

	return output;
}


ImageMat cvtColor(const ImageMat &origin, ImageMat::Type outputType)
{
	ImageMat output;
	if (origin.getType() == outputType) {
		output = origin;
		return output;
	}
	if (origin.getType() == ImageMat::Gray) {
		// Gray to non-gray is not supported
		throw std::runtime_error("Gray to non-gray is not supported");
	}
	if (outputType == ImageMat::Gray) {
		output = ImageMat(origin.getWidth(), origin.getHeight(), 1, outputType);
	}
	else {
		output = ImageMat(origin.getWidth(), origin.getHeight(), 3, outputType);
	}

	// Convert to BGR
	ImageMat bgrImage = cvtColorToBGR(origin);
	// Convert BGR To other
	uint32_t cols = bgrImage.getWidth();
	uint32_t rows = bgrImage.getHeight();
	ImageMat::Byte *originData = bgrImage.getRawData();
	ImageMat::Byte *outputData = output.getRawData();

	switch (outputType) {
	case ImageMat::Gray: {
		for (uint32_t i = 0; i < origin.getHeight(); ++i) {
			for (uint32_t j = 0; j < origin.getWidth(); ++j) {
				ImageMat::Byte *originPix = originData + i * cols * 3 + j * 3;
				outputData[i * cols + j] = static_cast<ImageMat::Byte>(originPix[2] * 0.299 + originPix[1] * 0.587 +
					originPix[0] * 0.114);
			}
		}
		break;
	}

	case ImageMat::YUV: {
		for (uint32_t i = 0; i < origin.getHeight(); ++i) {
			for (uint32_t j = 0; j < origin.getWidth(); ++j) {
				ImageMat::Byte *originPix = originData + i * cols * 3 + j * 3;
				ImageMat::Byte *outputPix = outputData + i * cols * 3 + j * 3;
				ImageMat::Byte R = originPix[2];
				ImageMat::Byte G = originPix[1];
				ImageMat::Byte B = originPix[0];
				// Y
				outputPix[0] = (ImageMat::Byte)(0.299 * R + 0.587 * G + 0.114 * B);
				// U
				outputPix[1] = (ImageMat::Byte)((-0.148 * R - 0.289 * G + 0.437 * B) / 0.437 * 0.5 + 128);
				// V
				outputPix[2] = (ImageMat::Byte)((0.615 * R - 0.515 * G - 0.100 * B) / 0.615 * 0.5 + 128);
			}
		}
		break;
	}

	case ImageMat::YCbCr: {
		for (uint32_t i = 0; i < origin.getHeight(); ++i) {
			for (uint32_t j = 0; j < origin.getWidth(); ++j) {
				ImageMat::Byte *originPix = originData + i * cols * 3 + j * 3;
				ImageMat::Byte *outputPix = outputData + i * cols * 3 + j * 3;
				ImageMat::Byte R = originPix[2];
				ImageMat::Byte G = originPix[1];
				ImageMat::Byte B = originPix[0];
				// Y
				outputPix[0] = (ImageMat::Byte)(0.299 * R + 0.587 * G + 0.114 * B);
				// Cb
				outputPix[1] = (ImageMat::Byte)(-0.1687 * R - 0.3313 * G + 0.500 * B + 128);
				// Cr
				outputPix[2] = (ImageMat::Byte)(0.5000 * R - 0.4187 * G - 0.0813 * B + 128);
			}
		}
		break;
	}

	case ImageMat::YIQ: {
		for (uint32_t i = 0; i < origin.getHeight(); ++i) {
			for (uint32_t j = 0; j < origin.getWidth(); ++j) {
				ImageMat::Byte *originPix = originData + i * cols * 3 + j * 3;
				ImageMat::Byte *outputPix = outputData + i * cols * 3 + j * 3;
				ImageMat::Byte R = originPix[2];
				ImageMat::Byte G = originPix[1];
				ImageMat::Byte B = originPix[0];
				// Y
				outputPix[0] = (ImageMat::Byte)(0.299 * R + 0.587 * G + 0.114 * B);
				// I
				outputPix[1] = (ImageMat::Byte)((0.596 * R - 0.274 * G - 0.322 * B) / 0.596 * 0.5 + 128);
				// Q
				outputPix[2] = (ImageMat::Byte)((0.211 * R - 0.523 * G + 0.312 * B) / 0.523 * 0.5 + 128);
			}
		}
		break;
	}

	case ImageMat::HSI: {
		for (uint32_t i = 0; i < origin.getHeight(); ++i) {
			for (uint32_t j = 0; j < origin.getWidth(); ++j) {
				ImageMat::Byte *originPix = originData + i * cols * 3 + j * 3;
				ImageMat::Byte *outputPix = outputData + i * cols * 3 + j * 3;
				ImageMat::Byte R = originPix[2];
				ImageMat::Byte G = originPix[1];
				ImageMat::Byte B = originPix[0];

				// H
				if (R == G && G == B) {
					outputPix[0] = 0;
				}
				else {
					double tmp = std::sqrt((R - G) * (R - G) + (R - B) * (G - B));
					double theta = (180 / M_PI) * acos(0.5 * (2 * R - G - B) / tmp);

					if (B <= G) {
						outputPix[0] = (ImageMat::Byte)(theta * 255.0 / 360);
					}
					else {
						outputPix[0] = (ImageMat::Byte)((360 - theta) * 255.0 / 360);
					}
				}
				
				// S
				ImageMat::Byte minSubpixel = R;
				if (G < minSubpixel) {
					minSubpixel = G;
				}
				if (B < minSubpixel) {
					minSubpixel = B;
				}
				ImageMat::Byte sum = R + G + B;

				if (sum != 0) {
					outputPix[1] = (ImageMat::Byte)((1 - 3.0 * minSubpixel / sum) * 255);
				}
				else {
					outputPix[1] = 0;
				}
				// I
				outputPix[2] = (ImageMat::Byte)(sum / 3);
			}
		}
		break;
	}

	case ImageMat::BGR:
		output = std::move(bgrImage);
		break;

	default:
		throw std::runtime_error("Unsupported Type");
		break;

	}
	output.setType(outputType);

	return output;
}