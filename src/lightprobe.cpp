/*
    This file is part of Tracer, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    [redacted] is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    [redacted] is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <tracer/lightprobe.h>
#include <tracer/vector.h>
#include <tracer/bitmap.h>

TRACER_NAMESPACE_BEGIN

static const int MultiplyDeBruijnBitPosition[64] =
{
	0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
	31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};

LightProbe::LightProbe(const std::string &filename): loaded(true) {
	Bitmap bitmap(filename);
	if (bitmap.cols() != bitmap.rows()) {
		throw TracerException("Width and height of a light probe must match!");
	}
	uint64_t size = bitmap.cols();
	if ((size & -size) != size) {
		throw TracerException("Size of a light probe must be power of 2!");
	}
	// Find the power of 2 in the size
	int M = MultiplyDeBruijnBitPosition[((uint32_t)(size * 0x077CB531U)) >> 27];

	// Convert to luminance and normalize
	Mipmap original = Mipmap(size, size);
	double sum = 0;
	for (int y = 0; y < size; ++y) {
		for (int x = 0; x < size; ++x) {
			double lum = bitmap.coeff(y, x).getLuminance();
			if (lum < 0) {
				lum = 0;
			}
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

TRACER_NAMESPACE_END
