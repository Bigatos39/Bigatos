#include "./HeaderFile/Geometry.h"
#include "./HeaderFile/Track.h"
#include "./HeaderFile/Simulator.h"

const float SigmaScatter = 0.01;
const float Sigma = 1;
const int det = 10;
const int d = 30;
const float Inf = 1000;

Simulator::Simulator(): fGeometry(), fScatter(SigmaScatter), fSigma(Sigma)
{
	fGeometry.zs.resize(det);
	for (int i = 0; i < det; i++)
	{
		fGeometry.zs[i] = d * i;
	}

	fRandom = new TRandom3();
}

void Simulator::Scatter(PositionParam &tp)
{
	const float alpha = fRandom -> Gaus(0, fScatter);
	tp.Tx() = tan(atan(tp.Tx()) + alpha);
}

void Simulator::Simulate(const PositionParam &param, Track &track)
{
	track.Points.clear();
	PositionParam curPoint = param;
	for (int i = 0; i < fGeometry.GetNDet(); i++)
	{
		const float z = fGeometry.zs[i];
		curPoint.X() = curPoint.X(z);
		curPoint.Z() = z;

		Scatter(curPoint);
		track.Points.push_back(curPoint);
	}

	track.hits.clear();
	for (int i = 0; i < fGeometry.GetNDet(); i++)
	{
		const float z = track.Points[i].z;
		float x = fRandom -> Gaus(track.Points[i].X(), fSigma);

		track.hits.push_back(Position(x, z));
	}
}

void Simulator::SetGeometry(const Geometry &g)
{
	fGeometry = g;
}

void Simulator::SetSigmaScatter(const float &s)
{
	fScatter = s;
}

void Simulator::SetSigma(const float &s)
{
	fSigma = s;
}
