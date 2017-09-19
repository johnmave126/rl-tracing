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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
	return Point2f(
		sample.x() < 0.5f ? (sqrtf(2 * sample.x()) - 1) : (1 - sqrtf(2 - 2 * sample.x())),
		sample.y() < 0.5f ? (sqrtf(2 * sample.y()) - 1) : (1 - sqrtf(2 - 2 * sample.y()))
	);
}

float Warp::squareToTentPdf(const Point2f &p) {
    return ((p.array() >= -1).all() && (p.array() <= 1).all()) ? (1 - abs(p.x())) * (1 - abs(p.y())) : 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
	Point2d polar = Point2d(
		sqrt(sample.y()),
		2 * M_PI * sample.x()
	);
	return Point2f(polar.x() * cos(polar.y()), polar.x() * sin(polar.y()));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
	return p.squaredNorm() < 1.0f - Epsilon ? 1.0 / M_PI : 0.0f;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
	Point2f euler = Point2f(asin(2 * sample.x() - 1), 2 * M_PI * sample.y());
	return Vector3f(
		cos(euler.x()) * cos(euler.y()),
		cos(euler.x()) * sin(euler.y()),
		sin(euler.x())
	);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
	return abs(v.squaredNorm() - 1.0f) < Epsilon ? 1.0 / (4 * M_PI) : 0.0f;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
	Point2f euler = Point2f(asin(sample.x()), 2 * M_PI * sample.y());
	return Vector3f(
		cos(euler.x()) * cos(euler.y()),
		cos(euler.x()) * sin(euler.y()),
		sin(euler.x())
	);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
	return abs(v.squaredNorm() - 1.0f) < Epsilon && v.z() >= 0 ? 1.0 / (2 * M_PI) : 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
	Point2f euler = Point2f(acos(1 - 2 * sample.x()) / 2, 2 * M_PI * sample.y());
	return Vector3f(
		sin(euler.x()) * cos(euler.y()),
		sin(euler.x()) * sin(euler.y()),
		cos(euler.x())
	);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
	return abs(v.squaredNorm() - 1.0f) < Epsilon && v.z() >= 0 ? v.dot(Vector3f(0, 0, 1)) / M_PI : 0.0f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
}

Point2f Warp::squareToLightProbe(const Point2f &sample, const LightProbe &probe) {
	Point2d local = sample.cast<double>();
	Point2d probed = Point2d(0.0f);
	int startCol = 0, startRow = 0;
	float baseUnit = 0.5f;
	for (int i = 0; i < probe.getCount(); i++, baseUnit /= 2.0f, startCol <<= 1, startRow <<= 1) {
		const LightProbe::Mipmap& map = probe.getMap(i);
		double horRatioUp = map.coeff(startRow, startCol) + map.coeff(startRow, startCol + 1);
		double horRatioBottom = map.coeff(startRow + 1, startCol) + map.coeff(startRow + 1, startCol + 1);
		double horRatio = horRatioUp / (horRatioUp + horRatioBottom);
		if (local.y() < horRatio) {
			local.y() /= horRatio;
			horRatio = horRatioUp;
		}
		else {
			probed.y() += baseUnit;
			local.y() = (local.y() - horRatio) / (1.0f - horRatio);
			horRatio = horRatioBottom;
			startRow++;
		}
		double vertRatio = map.coeff(startRow, startCol) / horRatio;
		if (local.x() < vertRatio) {
			local.x() /= vertRatio;
		}
		else {
			probed.x() += baseUnit;
			local.x() = (local.x() - vertRatio) / (1.0f - vertRatio);
			startCol++;
		}
	}
	return (probed + local * (2 * baseUnit)).cast<float>();
}

float Warp::squareToLightProbePdf(const Point2f & p, const LightProbe & probe)
{
	if (((p.array() < 0).any() || (p.array() >= 1).any())) {
		return 0.0f;
	}
	const LightProbe::Mipmap& topMap = probe.getMap(probe.getCount() - 1);
	return topMap.rows() * topMap.cols() * topMap.coeff(p.y() * topMap.rows(), p.x() * topMap.cols());
}

NORI_NAMESPACE_END
