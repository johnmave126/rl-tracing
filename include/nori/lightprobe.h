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

#pragma once

#include <nori/vector.h>

NORI_NAMESPACE_BEGIN

/// Mipmaps for a LightProbe image
class LightProbe {
public:
	typedef Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mipmap;
	/// Load an OpenEXR file with the specified filename
	LightProbe(const std::string &filename);

	/// Get the i-th mipmap of the probe
	const Mipmap& getMap(size_t idx) const {
		return mipmaps[idx];
	}

	/// Get the total number of maps
	size_t getCount() const {
		return mipmaps.size();
	}

private:
	std::vector<Mipmap> mipmaps;
};

NORI_NAMESPACE_END
