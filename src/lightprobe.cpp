/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/lightprobe.h>
#include <nori/vector.h>
#include <nori/bitmap.h>

NORI_NAMESPACE_BEGIN

static const int MultiplyDeBruijnBitPosition[64] =
{
	0,1,2,4,8,17,34,5,11,23,47,31,63,62,61,59,
	55,46,29,58,53,43,22,44,24,49,35,7,15,30,60,57,
	51,38,12,25,50,36,9,18,37,10,21,42,20,41,19,39,
	14,28,56,48,33,3,6,13,27,54,45,26,52,40,16,32
};

LightProbe::LightProbe(const std::string &filename) {
	Bitmap bitmap(filename);
	if (bitmap.cols() != bitmap.rows()) {
		throw NoriException("Width and height of a light probe must match!");
	}
	uint64_t size = bitmap.cols();
	if ((size & -size) != size) {
		throw NoriException("Size of a light probe must be power of 2!");
	}
	// Find the power of 2 in the size
	int M = MultiplyDeBruijnBitPosition[((uint64_t)(size * 0x022fdd63cc95386dULL)) >> 58];

	// Convert to luminance and normalize
	Mipmap original = Mipmap(size, size);
	float sum = 0;
	for (int y = 0; y < size; ++y) {
		for (int x = 0; x < size; ++x) {
			float lum = bitmap.coeff(y, x).getLuminance();
			original.coeffRef(y, x) = lum;
			sum += lum;
		}
	}
	original /= sum;
	mipmaps.push_back(std::move(original));

	for (int i = 1; i < M; i++) {
		size >>= 1;
		Mipmap downsample = Mipmap(size, size);
		const Mipmap& uppersample = mipmaps[i - 1];
		for (int y = 0; y < size; ++y) {
			for (int x = 0; x < size; ++x) {
				downsample.coeffRef(y, x) =
					uppersample.coeff(y << 1, x << 1)
					+ uppersample.coeff((y << 1), (x << 1) + 1)
					+ uppersample.coeff((y << 1) + 1, (x << 1))
					+ uppersample.coeff((y << 1) + 1, (x << 1) + 1);
			}
		}
		mipmaps.push_back(std::move(downsample));
	}
	std::reverse(std::begin(mipmaps), std::end(mipmaps));
}

NORI_NAMESPACE_END
