/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/Tools.h" // Has pi
#include <sstream> //std::ostringstream
#include <string.h> //strcasecmp
#include <cmath>

using namespace std;

namespace PLMD{
	namespace colvar{
		class QdistTwo : public Colvar {
			bool pbc;
			bool grid_mode_;
			bool serial;
			NeighborList *nl;
			bool invalidateList;
			bool firsttime;
			std::vector<double> active_q_;
			double cutoff;
			unsigned n_min_;
			double q_const_;
			std::vector<PLMD::AtomNumber> atomsToRequest;
			vector<AtomNumber> center_lista, around_species1_lista,around_species_lista;
			private:
			std::vector<Value*> valueDS;
			public:
			explicit QdistTwo(const ActionOptions&);
			~QdistTwo();
			// active methods:
			virtual void calculate();
			virtual void prepare();
			static void registerKeywords( Keywords& );

		};

		PLUMED_REGISTER_ACTION(QdistTwo,"QDISTTWO")

			void QdistTwo::registerKeywords( Keywords& keys ){
				Colvar::registerKeywords(keys);
				keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
				keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
				keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
				keys.add("optional","CUTOFF","cutoff distance, should not be bigger than half of the box size");
				keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
				keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
				keys.add("atoms","CENTER","Center atoms");
				keys.add("atoms","AROUND_SPECIES1","Around atoms 1");
				keys.add("optional","BOX_EDGE","set a reference value for the edge L of the simulation box. Compulsory for 'grid mode'");
				keys.add("optional","ACTIVE_Q","manually set which q frequencies will be considered");
				keys.add("optional","MAX","maximum grid value calculated. If LAMBDA is set, 2theta value is expected, otherwise q value. Default is reasonable when few primitive cells are simulated");
				keys.add("optional","MIN","minimum grid value calculated. If LAMBDA is set, 2theta value is expected, otherwise q value. Default is reasonable when few primitive cells are simulated");
				keys.add("optional","RESOLUTION","change q grid resolution, default is 3. Actual resolution depends on the size of the simulation box");
				keys.add("optional","NAME_PRECISION","set the number of digits used for components name");
				keys.addOutputComponent("ds","default","the instantaneous Debye Structure Factor at a given frequency q (or angle 2theta)");
				ActionWithValue::useCustomisableComponents(keys); //needed to have an unknown number of components
			}

		QdistTwo::QdistTwo(const ActionOptions&ao):
			PLUMED_COLVAR_INIT(ao),
			pbc(true),
			serial(false),
			invalidateList(true),
			firsttime(true)
		{
			parseFlag("SERIAL",serial);
			parse("CUTOFF",cutoff);


			bool nopbc=!pbc;
			parseFlag("NOPBC",nopbc);
			pbc=!nopbc;
			// pair stuff
			bool dopair=false;
			parseFlag("PAIR",dopair);

			// neighbor list stuff
			bool doneigh=false;
			double nl_cut=0.0;
			int nl_st=0;
			parseFlag("NLIST",doneigh);
			if(doneigh){
				parse("NL_CUTOFF",nl_cut);
				if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
				parse("NL_STRIDE",nl_st);
				if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
			}
			parseVector("ACTIVE_Q",active_q_);

			if (active_q_.size()>0)
			{
				//print info
				log.printf("  using only %d manually selected q:",active_q_.size());
				for (unsigned q=0; q<active_q_.size(); q++)
					log.printf("  %g",active_q_[q]);
				log.printf("\n");
				grid_mode_=false;
			}
			else
				grid_mode_=true;
			//////- get q frequency grid
			double box_edge=-1;
			double max=-1;
			double min=-1;
			unsigned resolution=0;
			parse("BOX_EDGE",box_edge);
			parse("MIN",min);
			parse("MAX",max);
			parse("RESOLUTION",resolution);
			if (grid_mode_)
			{
				if (resolution==0)
					resolution=3; //default resolution
				//set the grid
				q_const_=PLMD::pi/box_edge/resolution;
				if (min==-1)
					min=4*PLMD::pi/box_edge; //below this is not physical
				n_min_=std::ceil(min/q_const_);
				unsigned n_max;
				if (max==-1)
					n_max=n_min_+149; //good guess for small simulations
				else
				{
					n_max=std::floor(max/q_const_);
				}
				active_q_.resize(n_max-n_min_+1,q_const_);
				for (unsigned q=0; q<active_q_.size(); q++)
					active_q_[q]*=(q+n_min_);
				//print grid info
				log.printf("  using a grid on q space:\n");
				log.printf("    reference BOX_EDGE L = %g\n",box_edge);
				log.printf("    MIN q = %g [should be greater than 2pi/L=%g]\n",q_const_*n_min_,2*PLMD::pi/box_edge);
				log.printf("    MAX q = %g\n",q_const_*n_max);
				log.printf("    q RESOLUTION = %d --> %d grid points\n",resolution,active_q_.size()); 
			}
				else
					plumed_massert(box_edge==-1 && min==-1 && max==-1 && resolution==0,"if specific ACTIVE_Q are given, no grid parameter (BOX_EDGE,MIN,MAX,RESOLUTION) should be set");
			//////add colvar components
			valueDS.resize(active_q_.size());
			std::string val="q";
			std::ostringstream oss;
			unsigned name_precision=7;
			parse("NAME_PRECISION",name_precision);
			oss.precision(name_precision);
			log.printf("  components name are %s values, with NAME_PRECISION = %d\n",val.c_str(),name_precision);
			for (unsigned q=0; q<active_q_.size(); q++)
			{ 
				oss.str("");
				oss<<"ds-"<<active_q_[q];
                                addComponentWithDerivatives(oss.str());
                                componentIsNotPeriodic(oss.str());
				valueDS[q]=getPntrToComponent(oss.str());
			}
			//////////
                        parseAtomList("CENTER",center_lista);
			parseAtomList("AROUND_SPECIES1",around_species1_lista);
			around_species_lista.reserve ( around_species1_lista.size() );
			around_species_lista.insert (  around_species_lista.end() , around_species1_lista.begin(),  around_species1_lista.end() );
			if(center_lista.size()>0 && around_species1_lista.size()>0){
				if(doneigh)  nl= new NeighborList(center_lista,around_species_lista,serial,dopair,pbc,getPbc(),comm,nl_cut,nl_st);
				else         nl= new NeighborList(center_lista,around_species_lista,serial,dopair,pbc,getPbc(),comm);
			} else {
				error("CENTER, AROUND_SPECIES1 should be explicitly defined.");
			}
			atomsToRequest.reserve (center_lista.size() + around_species_lista.size() );
			atomsToRequest.insert (atomsToRequest.end(), center_lista.begin(), center_lista.end() );
			atomsToRequest.insert (atomsToRequest.end(), around_species_lista.begin(), around_species_lista.end() );
			requestAtoms(atomsToRequest);

			log.printf("  first group:\n");
			for(unsigned int i=0;i<center_lista.size();++i){
				if ( (i+1) % 25 == 0 ) log.printf("  \n");
				log.printf("  %d", center_lista[i].serial());
			}
			log.printf("  \n  second group:\n");
			for(unsigned int i=0;i<around_species1_lista.size();++i){
				if ( (i+1) % 25 == 0 ) log.printf("  \n");
				log.printf("  %d", around_species1_lista[i].serial());
			}
			log.printf("  \n");
			if(pbc) log.printf("  using periodic boundary conditions\n");
			else    log.printf("  without periodic boundary conditions\n");
			

                          checkRead(); 
		}

		QdistTwo::~QdistTwo(){
			delete nl;
		}

		void QdistTwo::prepare(){
			if(nl->getStride()>0){
				requestAtoms(atomsToRequest);
				if(firsttime || (getStep()%nl->getStride()==0)){
					invalidateList=true;
					firsttime=false;
				}else{
					invalidateList=false;
					if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
				}
				if(getExchangeStep()) firsttime=true;
			}
		}

		// calculator
		void QdistTwo::calculate()
		{
			if(nl->getStride()>0 && invalidateList){
				nl->update(getPositions());
			}
			// Print positions of A atoms
			unsigned stride=comm.Get_size();
			unsigned rank=comm.Get_rank();
			if(serial){
				stride=1;
				rank=0;
			}else{
				stride=comm.Get_size();
				rank=comm.Get_rank();
			}
			const unsigned nn=center_lista.size();
			std::vector<double> DebyeS(active_q_.size(),0.);
			// Loop over center atoms
			for(unsigned int i=rank;i<nn;i+=stride) {
				std::vector<double> DebyeS_group(active_q_.size(),0.);
				// Loop over around species
				std::vector<unsigned> neighbors;
				neighbors=nl->getNeighbors(i);
				for(unsigned int j=0;j<neighbors.size();j+=1) {   
					Vector dij; 
					double window=1;
					if (atomsToRequest[i].serial()!=atomsToRequest[neighbors[j]].serial()) { 
						// i different from j
						dij = pbcDistance(getPosition(i),getPosition(neighbors[j]));
       if ( (neighbors[j]>=center_lista.size()) && (neighbors[j]<(center_lista.size()+around_species1_lista.size())) ) {
						double R_ij = dij.modulo();
						if (R_ij>=cutoff) 
							continue;
						const double window_arg=PLMD::pi*R_ij/cutoff;
						window=sin(window_arg)/window_arg;
						double base_cos; //this throws a warning, but is correct (and faster)
						double base_sin;
						double prev_cos;
						double prev_sin;
						base_cos=cos(q_const_*R_ij);
						base_sin=sin(q_const_*R_ij);
						prev_cos=cos(q_const_*R_ij*(n_min_-1));//room for improvements?
						prev_sin=sin(q_const_*R_ij*(n_min_-1));
						for (unsigned q=0; q<active_q_.size(); q++)
						{
							const double QR=active_q_[q]*R_ij;
							double cosQR;
							double sinQR;
							cosQR=base_cos*prev_cos-base_sin*prev_sin;
							sinQR=base_cos*prev_sin+base_sin*prev_cos;
							prev_cos=cosQR;
							prev_sin=sinQR;
							DebyeS_group[q]+=window*sinQR/QR;
						}

					}
				}
}
				for (unsigned q=0; q<active_q_.size(); q++)
					DebyeS[q]+=DebyeS_group[q]/nn;
			}
			if(!serial){
				comm.Sum(DebyeS);
			}
			for (unsigned q=0; q<active_q_.size(); q++)
			{
				valueDS[q]->set(1+DebyeS[q]);
			}

		}
	}
}	
