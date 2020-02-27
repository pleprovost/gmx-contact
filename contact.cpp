/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include <string>
#include <vector>
#include <future>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <functional>
#include <gromacs/trajectoryanalysis.h>

using namespace gmx;

/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */
class AnalysisTemplate : public TrajectoryAnalysisModule
{
public:
    AnalysisTemplate();

    virtual void initOptions(IOptionsContainer          *options,
			     TrajectoryAnalysisSettings *settings);
    virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
			      const TopologyInformation        &top);

    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
			      TrajectoryAnalysisModuleData *pdata);

    virtual void finishAnalysis(int nframes);
    virtual void writeOutput();

private:
    class ModuleData;

    std::string                      fnDist_;
    double                           cutoff_;
    Selection                        refsel_;
    Selection                        sel_;

    AnalysisData                     data_;
    AnalysisDataAverageModulePointer avem_;
    
    std::vector<int> contactMatrix_;
    std::ofstream outputStream_;
    size_t refSize_;
    size_t selSize_;

    int nb_total_contact_ = 0;
};


AnalysisTemplate::AnalysisTemplate()
    : cutoff_(1.0)
{
    registerAnalysisDataset(&data_, "avedist");
    cutoff_ *= cutoff_;
}


void
AnalysisTemplate::initOptions(IOptionsContainer          *options,
                              TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "This is a template for writing your own analysis tools for",
        "GROMACS. The advantage of using GROMACS for this is that you",
        "have access to all information in the topology, and your",
        "program will be able to handle all types of coordinates and",
        "trajectory files supported by GROMACS. In addition,",
        "you get a lot of functionality for free from the trajectory",
        "analysis library, including support for flexible dynamic",
        "selections. Go ahead an try it![PAR]",
        "To get started with implementing your own analysis program,",
        "follow the instructions in the README file provided.",
        "This template implements a simple analysis programs that calculates",
        "average distances from a reference group to one or more",
        "analysis groups."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("o")
		       .filetype(eftPlot).outputFile()
		       .store(&fnDist_).defaultBasename("avedist")
		       .description("Average distances from reference group"));

    options->addOption(SelectionOption("reference")
		       .store(&refsel_).required()
		       .description("Reference group to calculate distances from"));
    options->addOption(SelectionOption("select")
		       .store(&sel_).required()
		       .description("Group to calculate distances to"));

    options->addOption(DoubleOption("cutoff").store(&cutoff_)
		       .description("Cutoff for distance calculation (0 = no cutoff)"));

    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}


void
AnalysisTemplate::initAnalysis(const TrajectoryAnalysisSettings &settings,
                               const TopologyInformation         & /*top*/)
{
    data_.setColumnCount(0, 1); 
    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);

    refSize_ = refsel_.posCount();
    selSize_ = sel_.posCount();

    contactMatrix_.resize(refSize_*selSize_);
    std::fill(contactMatrix_.begin(), contactMatrix_.end(), 0);
    
    if (!fnDist_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
	    new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnDist_);
        plotm->setTitle("Average distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        data_.addModule(plotm);
    }
}

void contact(int start, int stop, const Selection &ref, const Selection &sel,
	     std::vector<int> contactMatrix, int& nb_contact)
{

    
}

void
AnalysisTemplate::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle         dh     = pdata->dataHandle(data_);
    const Selection           &refsel = pdata->parallelSelection(refsel_);
    const Selection           &sel = pdata->parallelSelection(sel_);
    
    dh.startFrame(frnr, fr.time);
    int nb_contact = 0;

    auto l = [&] (int start, int stop) -> int {
	int count=0;
	for (size_t i = start; i < stop; ++i)
	{
	    for (size_t j = 0; j < selSize_; ++j)
	    {
		SelectionPosition p = refsel.position(i);
		SelectionPosition q = sel.position(j);
		if ( distance2(p.x(), q.x()) < cutoff_)
		{
		    contactMatrix_.at(refSize_*i+j) += 1;
		    count++;
		}
	    }
	}
	return count;
    };

    std::vector<std::future<int>> futures{};
    int n = std::thread::hardware_concurrency();
    int chunksize = refSize_ / n;
    
    std::vector<std::pair<int, int>> startstop;
    for ( int k = 0; k < n; k++ )
    {
	startstop.push_back(std::pair<int, int>(k*chunksize, (k+1)*chunksize));
    }
    
    // Launch Futures
    for ( auto& pair : startstop )
    {
	futures.push_back(std::async(std::launch::async, l,
				     pair.first,
				     pair.second));
    }
    // Get result
    for ( auto &fut : futures )
    {
	nb_contact += fut.get();
    }

    // for (size_t i = 0; i < refSize_; ++i)
    // {
    //     for (size_t j = 0; j < selSize_; ++j)
    //     {
    // 	SelectionPosition p = refsel.position(i);
    // 	SelectionPosition q = sel.position(j);
    // 	if ( distance2(p.x(), q.x()) < cutoff_)
    // 	{
    // 	    contactMatrix_.at(refSize_*i+j) += 1;
    // 	    nb_contact++;
    // 	}
    //     }
    // }
	
    dh.setPoint(0, nb_contact);
    nb_total_contact_ += nb_contact;
    dh.finishFrame();
}


void
AnalysisTemplate::finishAnalysis(int nframes)
{
    
}


void
AnalysisTemplate::writeOutput()
{
    // Write tjhe output contact matrix
    outputStream_.open("contact_matrix.dat");
    for (size_t i = 0; i < refSize_; ++i)
    { 
   	for (size_t j = 0; j < selSize_; ++j)
   	{
   	    outputStream_ << std::setw(6)
   			  << 1.0*contactMatrix_.at(refSize_*i+j)/nb_total_contact_ << " "; 
   	}
   	outputStream_ << "\n";
    }
    outputStream_.flush();
    outputStream_.close();
}

/*! \brief
 * The main function for the analysis template.
 */
int
main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<AnalysisTemplate>(argc, argv);
}
