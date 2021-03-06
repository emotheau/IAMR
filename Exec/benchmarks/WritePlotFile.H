#ifndef _WritePlotFile_H_
#define _WritePlotFile_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_AmrData.H>

// Write a MultiFab to a plotfile named "name" using the HyperCLawV1.1 parallel
// plotfile data format
void
writePlotFile (const char*     name,
	       const MultiFab& mf,
	       const Geometry& geom,
	       const IntVect&  refRatio,
	       Real            bgVal = 0.0);

// Write an array of multifabs pulling all the data from the AmrData
void
WritePlotFile(const Vector<MultiFab*> mfa,
	      AmrData&               amrdToMimic,
	      const std::string&     oFile,
	      bool                   verbose);


#endif

