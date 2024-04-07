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
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/Tools.h" // Has pi
#include "tools/Angle.h"
#include "tools/SwitchingFunction.h"
#include "tools/HistogramBead.h"
#include <cmath>

using namespace std;

namespace PLMD{
	namespace colvar{
		class dsfAll : public Colvar {
			bool pbc;
			bool serial;
			NeighborList *nl;
			bool invalidateList;
			bool firsttime;
			double active_q1;
			double active_q2;
			double active_q3;
			double cutoff1;
			double cutoff2;
			double cutoff3;
			std::vector<PLMD::AtomNumber> atomsToRequest;
			SwitchingFunction switchingFunction1;
			SwitchingFunction switchingFunction2;
			SwitchingFunction switchingFunction3;
			SwitchingFunction distSwitchingFunction1;
			SwitchingFunction distSwitchingFunction2;
			SwitchingFunction distSwitchingFunction3;
			SwitchingFunction coordSwitchingFunction1;
			SwitchingFunction coordSwitchingFunction2;
			SwitchingFunction coordSwitchingFunction3;
			vector<AtomNumber> center_lista, around_species1_lista, around_species2_lista, around_species3_lista, around_species_lista ;
			public:
			explicit dsfAll(const ActionOptions&);
			~dsfAll();
			// active methods:
			virtual void calculate();
			virtual void prepare();
			static void registerKeywords( Keywords& keys );

		};

		PLUMED_REGISTER_ACTION(dsfAll,"DSFALL")

			void dsfAll::registerKeywords( Keywords& keys ){
				Colvar::registerKeywords(keys);
				keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
				keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
				keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
				keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
				keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
				keys.add("atoms","CENTER","Center atoms");
				keys.add("atoms","AROUND_SPECIES1","Around atoms 1");
				keys.add("atoms","AROUND_SPECIES2","Around atoms 2");
				keys.add("atoms","AROUND_SPECIES3","Around atoms 3");
				keys.add("optional","SWITCH1","Switching function 1.");
				keys.add("optional","SWITCH2","Switching function 2.");
				keys.add("optional","SWITCH3","Switching function 3.");
				keys.add("optional","DIST_SWITCH1","Switching function 1.");
				keys.add("optional","DIST_SWITCH2","Switching function 2.");
				keys.add("optional","DIST_SWITCH3","Switching function 3.");
				keys.add("optional","COORD_SWITCH1","Coordination switching function 1.");
				keys.add("optional","COORD_SWITCH2","Coordination switching function 2.");
				keys.add("optional","COORD_SWITCH3","Coordination switching function 3.");
				keys.add("optional","ACTIVE_Q1","manually set which q frequencies will be considered");
				keys.add("optional","ACTIVE_Q2","manually set which q frequencies will be considered");
				keys.add("optional","ACTIVE_Q3","manually set which q frequencies will be considered");
				keys.add("optional","CUTOFF1","cutoff distance, should not be bigger than half of the box size");
				keys.add("optional","CUTOFF2","cutoff distance, should not be bigger than half of the box size");
				keys.add("optional","CUTOFF3","cutoff distance, should not be bigger than half of the box size");
				keys.add("optional","BOX_EDGE","set a reference value for the edge L of the simulation box. Compulsory for 'grid mode'");
				keys.add("optional","MAX","maximum grid value calculated. If LAMBDA is set, 2theta value is expected, otherwise q value. Default is reasonable when few primitive cells are simulated");
				keys.add("optional","MIN","minimum grid value calculated. If LAMBDA is set, 2theta value is expected, otherwise q value. Default is reasonable when few primitive cells are simulated");
				keys.add("optional","RESOLUTION","change q grid resolution, default is 3. Actual resolution depends on the size of the simulation box");
			}

		dsfAll::dsfAll(const ActionOptions&ao):
			PLUMED_COLVAR_INIT(ao),
			pbc(true),
			serial(false),
			invalidateList(true),
			firsttime(true)
		{
			parseFlag("SERIAL",serial);
			parseAtomList("CENTER",center_lista);
			parseAtomList("AROUND_SPECIES1",around_species1_lista);
			parseAtomList("AROUND_SPECIES2",around_species2_lista);
			parseAtomList("AROUND_SPECIES3",around_species3_lista);
			around_species_lista.reserve ( around_species1_lista.size() + around_species2_lista.size() +around_species3_lista.size() );
			around_species_lista.insert (  around_species_lista.end() , around_species1_lista.begin(),  around_species1_lista.end() );
			around_species_lista.insert (  around_species_lista.end() , around_species2_lista.begin(),  around_species2_lista.end() );
			around_species_lista.insert (  around_species_lista.end() , around_species3_lista.begin(),  around_species3_lista.end() );
			parse("ACTIVE_Q1",active_q1);
			parse("ACTIVE_Q2",active_q2);
			parse("ACTIVE_Q3",active_q3);
			parse("CUTOFF1",cutoff1);
			parse("CUTOFF2",cutoff2);
			parse("CUTOFF3",cutoff3);


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
			addValueWithDerivatives(); setNotPeriodic();
			if(center_lista.size()>0 && around_species1_lista.size()>0  && around_species2_lista.size()>0 && around_species3_lista.size()>0){
				if(doneigh)  nl= new NeighborList(center_lista,around_species_lista,serial,dopair,pbc,getPbc(),comm,nl_cut,nl_st);
				else         nl= new NeighborList(center_lista,around_species_lista,serial,dopair,pbc,getPbc(),comm);
			} else {
				error("CENTER, AROUND_SPECIES1, AROUND_SPECIES2, and AROUND_SPECIES3 should be explicitly defined.");
			}
			atomsToRequest.reserve ( center_lista.size() + around_species_lista.size() );
			atomsToRequest.insert (atomsToRequest.end(), center_lista.begin(), center_lista.end() );
			atomsToRequest.insert (atomsToRequest.end(), around_species_lista.begin(), around_species_lista.end() );

			requestAtoms(atomsToRequest);
			log.printf("  between four groups of %u, %u, %u and %u atoms\n",static_cast<unsigned>(center_lista.size()),static_cast<unsigned>(around_species1_lista.size()),static_cast<unsigned>(around_species2_lista.size()), static_cast<unsigned>(around_species3_lista.size()) );
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
			log.printf("  \n  third group:\n");
			for(unsigned int i=0;i<around_species2_lista.size();++i){
				if ( (i+1) % 25 == 0 ) log.printf("  \n");
				log.printf("  %d", around_species2_lista[i].serial());
			}
			log.printf("  \n  fourth group:\n");
			for(unsigned int i=0;i<around_species3_lista.size();++i){
				if ( (i+1) % 25 == 0 ) log.printf("  \n");
				log.printf("  %d", around_species3_lista[i].serial());
			}
			log.printf("  \n");
			if(pbc) log.printf("  using periodic boundary conditions\n");
			else    log.printf("  without periodic boundary conditions\n");
			if(dopair) log.printf("  with PAIR option\n");
			if(doneigh){
				log.printf("  using neighbor lists with\n");
				log.printf("  update every %d steps and cutoff %f \n",nl_st,nl_cut);
			}
			string sw,errors;

			// Intensity switching functions
			parse("SWITCH1",sw);
			if(sw.length()>0){
				switchingFunction1.set(sw,errors);
				if( errors.length()!=0 ) error("problem reading SWITCH1 keyword : " + errors );
			} else {
				error("No switching function 1 defined. Please define SWITCH1!");
			}
			parse("SWITCH2",sw);
			if(sw.length()>0){
				switchingFunction2.set(sw,errors);
				if( errors.length()!=0 ) error("problem reading SWITCH2 keyword : " + errors );
			} else {
				error("No switching function 2 defined. Please define SWITCH2!");
			}
			parse("SWITCH3",sw);
			if(sw.length()>0){
				switchingFunction3.set(sw,errors);
				if( errors.length()!=0 ) error("problem reading SWITCH3 keyword : " + errors );
			} else {
				error("No switching function 3 defined. Please define SWITCH3!");
			}
			// Distance switching functions
			parse("DIST_SWITCH1",sw);
			if(sw.length()>0){
				distSwitchingFunction1.set(sw,errors);
				if( errors.length()!=0 ) error("problem reading COORD_SWITCH1 keyword : " + errors );
			} else {
				error("No coordination switching function 1 defined. Please define COORD_SWITCH1!");
			}
			parse("DIST_SWITCH2",sw);
			if(sw.length()>0){
				distSwitchingFunction2.set(sw,errors);
				if( errors.length()!=0 ) error("problem reading COORD_SWITCH2 keyword : " + errors );
			} else {
				error("No coordination switching function 2 defined. Please define COORD_SWITCH2!");
			}
			parse("DIST_SWITCH3",sw);
			if(sw.length()>0){
				distSwitchingFunction3.set(sw,errors);
				if( errors.length()!=0 ) error("problem reading COORD_SWITCH3 keyword : " + errors );
			} else {
				error("No coordination switching function 3 defined. Please define COORD_SWITCH3!");
			}
			// Coordination switching functions
			parse("COORD_SWITCH1",sw);
			if(sw.length()>0){
				coordSwitchingFunction1.set(sw,errors);
				if( errors.length()!=0 ) error("problem reading COORD_SWITCH1 keyword : " + errors );
			} else {
				error("No coordination switching function 1 defined. Please define COORD_SWITCH1!");
			}
			parse("COORD_SWITCH2",sw);
			if(sw.length()>0){
				coordSwitchingFunction2.set(sw,errors);
				if( errors.length()!=0 ) error("problem reading COORD_SWITCH2 keyword : " + errors );
			} else {
				error("No coordination switching function 2 defined. Please define COORD_SWITCH2!");
			}
			parse("COORD_SWITCH3",sw);
			if(sw.length()>0){
				coordSwitchingFunction3.set(sw,errors);
				if( errors.length()!=0 ) error("problem reading COORD_SWITCH3 keyword : " + errors );
			} else {
				error("No coordination switching function 3 defined. Please define COORD_SWITCH3!");
			}

			checkRead();
		}

		dsfAll::~dsfAll(){
			delete nl;
		}

		void dsfAll::prepare(){
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
		void dsfAll::calculate()
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
			double CVvalue=0.;
			vector<Vector> deriv(getNumberOfAtoms());
			Tensor virial;
			// Loop over center atoms
			for(unsigned int i=rank;i<nn;i+=stride) {
				// Loop over around species
				double DS1=0., DS2=0., DS3=0.;
				vector<Vector> deriv1s(getNumberOfAtoms()), deriv2s(getNumberOfAtoms()), deriv3s(getNumberOfAtoms());
				Tensor virial1s, virial2s, virial3s;
				double coordination1=0., coordination2=0., coordination3=0.;
				vector<Vector> deriv1c(getNumberOfAtoms()), deriv2c(getNumberOfAtoms()), deriv3c(getNumberOfAtoms());
				Tensor virial1c, virial2c, virial3c;
				std::vector<unsigned> neighbors;
				neighbors=nl->getNeighbors(i);
				for(unsigned int j=0;j<neighbors.size();j+=1) {   
					Vector dij;
					double dfunc; 
					if (atomsToRequest[i].serial()!=atomsToRequest[neighbors[j]].serial()) { 
						// i different from j
						dij = pbcDistance(getPosition(i),getPosition(neighbors[j]));
						double R_ij = dij.modulo();
						if ( (neighbors[j]>=center_lista.size()) && (neighbors[j]<(center_lista.size()+around_species1_lista.size())) ) {
							if (R_ij>=cutoff1) 
								continue;
							double window=1.;
							double d_window=0.;
							const double window_arg=PLMD::pi*R_ij/cutoff1;
							window=std::sin(window_arg)/window_arg;
							d_window=(window_arg*std::cos(window_arg)-std::sin(window_arg))/window_arg/R_ij;
							const double QR1=active_q1*R_ij;
							DS1+=std::sin(QR1)/QR1*window;
							Vector dds = dij*((window*(QR1*std::cos(QR1)-std::sin(QR1))/R_ij+d_window*std::sin(QR1))/QR1/R_ij);
							deriv1s[i] -= dds;
							deriv1s[neighbors[j]] += dds;  
							Tensor vvs(dij,dds);
							virial1s -= vvs;
							// Atom belongs to around species 1
							coordination1 += distSwitchingFunction1.calculateSqr(dij.modulo2(),dfunc);
							Vector ddc(dfunc*dij);
							Tensor vvc(ddc,dij);
							virial1c -= vvc;
							deriv1c[i]-=ddc;
							deriv1c[neighbors[j]]+=ddc;

						} else if ( (neighbors[j]>=(center_lista.size()+around_species1_lista.size()) ) && (neighbors[j]<(center_lista.size()+around_species1_lista.size()+around_species2_lista.size())) ){
							if (R_ij>=cutoff2) 
								continue;
							double window=1.;
							double d_window=0.;
							const double window_arg=PLMD::pi*R_ij/cutoff2;
							window=std::sin(window_arg)/window_arg;
							d_window=(window_arg*std::cos(window_arg)-std::sin(window_arg))/window_arg/R_ij;
							const double QR2=active_q2*R_ij;
							DS2+=std::sin(QR2)/QR2*window;
							Vector dds = dij*((window*(QR2*std::cos(QR2)-std::sin(QR2))/R_ij+d_window*std::sin(QR2))/QR2/R_ij);
							deriv2s[i] -= dds;
							deriv2s[neighbors[j]] += dds;  
							Tensor vvs(dij,dds);
							virial2s -= vvs;
							// Atom belongs to around species 2
							coordination2 += distSwitchingFunction2.calculateSqr(dij.modulo2(),dfunc);
							Vector ddc(dfunc*dij);
							Tensor vvc(ddc,dij);
							virial2c -= vvc;
							deriv2c[i]-=ddc;
							deriv2c[neighbors[j]]+=ddc;

						} else if ( neighbors[j]>=(center_lista.size()+around_species1_lista.size()+around_species2_lista.size()) ){
							if (R_ij>=cutoff3) 
								continue;
							double window=1.;
							double d_window=0.;
							const double window_arg=PLMD::pi*R_ij/cutoff3;
							window=std::sin(window_arg)/window_arg;
							d_window=(window_arg*std::cos(window_arg)-std::sin(window_arg))/window_arg/R_ij;
							const double QR3=active_q3*R_ij;
							DS3+=std::sin(QR3)/QR3*window;
							Vector dds = dij*((window*(QR3*std::cos(QR3)-std::sin(QR3))/R_ij+d_window*std::sin(QR3))/QR3/R_ij);
							deriv3s[i] -= dds;
							deriv3s[neighbors[j]] += dds;  
							Tensor vvs(dij,dds);
							virial3s -= vvs;
							// Atom belongs to around species 3
							coordination3 += distSwitchingFunction3.calculateSqr(dij.modulo2(),dfunc);
							Vector ddc(dfunc*dij);
							Tensor vvc(ddc,dij);
							virial3c -= vvc;
							deriv3c[i]-=ddc;
							deriv3c[neighbors[j]]+=ddc;

						}			
					}
				}
				// CVvalue+=DS/nn;
				double dfunc1s,dfunc2s,dfunc3s;
				DS1 = DS1+1;
				DS2 = DS2+1;
				DS3 = DS3+1;
				double CVvalue1s = switchingFunction1.calculate(DS1,dfunc1s);
				double CVvalue2s = switchingFunction2.calculate(DS2,dfunc2s);
				double CVvalue3s = switchingFunction3.calculate(DS3,dfunc3s);
				dfunc1s *= DS1;
				dfunc2s *= DS2;
				dfunc3s *= DS3;
				double dfunc1c,dfunc2c,dfunc3c;
				double CVvalue1c = coordSwitchingFunction1.calculate(coordination1,dfunc1c);
				double CVvalue2c = coordSwitchingFunction2.calculate(coordination2,dfunc2c);
				double CVvalue3c = coordSwitchingFunction3.calculate(coordination3,dfunc3c);
				dfunc1c *= coordination1;
				dfunc2c *= coordination2;
				dfunc3c *= coordination3;

				CVvalue += CVvalue1s*CVvalue2s*CVvalue3s*CVvalue1c*CVvalue2c*CVvalue3c;
				for(unsigned int j=0;j<getNumberOfAtoms();j+=1) {   
					for(unsigned int k=0;k<3;k+=1) {  
						deriv[j][k] += deriv1s[j][k]*dfunc1s*CVvalue2s*CVvalue3s*CVvalue1c*CVvalue2c*CVvalue3c + deriv2s[j][k]*dfunc2s*CVvalue1s*CVvalue3s*CVvalue1c*CVvalue2c*CVvalue3c + deriv3s[j][k]*dfunc3s*CVvalue2s*CVvalue1s*CVvalue1c*CVvalue2c*CVvalue3c + deriv1c[j][k]*dfunc1c*CVvalue2c*CVvalue3c*CVvalue1s*CVvalue2s*CVvalue3s +  deriv2c[j][k]*dfunc2c*CVvalue1c*CVvalue3c*CVvalue1s*CVvalue2s*CVvalue3s + deriv3c[j][k]*dfunc3c*CVvalue1c*CVvalue2c*CVvalue1s*CVvalue2s*CVvalue3s;
					}
				}
				for(unsigned int j=0;j<3;j+=1) {   
					for(unsigned int k=0;k<3;k+=1) {
						virial[j][k] += virial1s[j][k]*dfunc1s*CVvalue2s*CVvalue3s*CVvalue1c*CVvalue2c*CVvalue3c + virial2s[j][k]*dfunc2s*CVvalue1s*CVvalue3s*CVvalue1c*CVvalue2c*CVvalue3c + virial3s[j][k]*dfunc3s*CVvalue2s*CVvalue1s*CVvalue1c*CVvalue2c*CVvalue3c + virial1c[j][k]*dfunc1c*CVvalue2c*CVvalue3c*CVvalue1s*CVvalue2s*CVvalue3s +  virial2c[j][k]*dfunc2c*CVvalue1c*CVvalue3c*CVvalue1s*CVvalue2s*CVvalue3s + virial3c[j][k]*dfunc3c*CVvalue1c*CVvalue2c*CVvalue1s*CVvalue2s*CVvalue3s;
					}
				}
			}
			if(!serial){
				comm.Sum(CVvalue);
				if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());
				comm.Sum(virial);
			}
			setValue(CVvalue);
			for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]);
			setBoxDerivatives  (virial);

		}
	}
}	
