/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/Tools.h" // Has pi
#include "tools/SwitchingFunction.h"
#include <string>

using namespace std;

namespace PLMD{
namespace colvar{

class Oct : public Colvar {
  bool pbc;
  bool serial;
  NeighborList *nl;
  bool invalidateList;
  bool firsttime;
  std::vector<PLMD::AtomNumber> atomsToRequest;
  SwitchingFunction switchingFunction1;
  SwitchingFunction switchingFunction2;
  SwitchingFunction switchingFunction3;
  SwitchingFunction coordSwitchingFunction1;
  SwitchingFunction coordSwitchingFunction2;
  SwitchingFunction coordSwitchingFunction3;
  vector<AtomNumber> center_lista, around_species1_lista, around_species2_lista, around_species3_lista, around_species_lista;
  mutable PLMD::OFile outputFile;
public:
  explicit Oct(const ActionOptions&);
  ~Oct();
// active methods:
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords( Keywords& keys );
  double dotprod(Vector,Vector);
  double norm2(Vector);
};

PLUMED_REGISTER_ACTION(Oct,"OCT")

void Oct::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list 1");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","CENTER","Center atoms");
  keys.add("atoms","AROUND_SPECIES1","Around atoms 1");
  keys.add("atoms","AROUND_SPECIES2","Around atoms 1");
  keys.add("atoms","AROUND_SPECIES3","Around atoms 1");
  keys.add("optional","SWITCH1","Switching function 1.");
  keys.add("optional","SWITCH2","Switching function 1.");
  keys.add("optional","SWITCH3","Switching function 1.");
  keys.add("optional","COORD_SWITCH1","Coordination switching function 1.");
  keys.add("optional","COORD_SWITCH2","Coordination switching function 1.");
  keys.add("optional","COORD_SWITCH3","Coordination switching function 1.");
}

Oct::Oct(const ActionOptions&ao):
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
  // Distance switching functions
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

 outputFile.link(*this);
 outputFile.open("output.xyz");
  checkRead();
}

Oct::~Oct(){
  delete nl;
}

void Oct::prepare(){
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
//OFile outputFile;
// calculator
void Oct::calculate()
{
 if(nl->getStride()>0 && invalidateList){
   nl->update(getPositions());
 }
 unsigned stride=comm.Get_size();
 unsigned rank=comm.Get_rank();
 if(serial){
   stride=1;
   rank=0;
 }else{
   stride=comm.Get_size();
   rank=comm.Get_rank();
 }

//OFile outputFile;
//outputFile.link(*this);
double CVvalue_edge=0.;
double CVvalue_face=0.;
double CVvalue_corner=0.;
//std::map<int, int> coordi_map;
//std::map<int,int>::iterator it;
const unsigned nn=center_lista.size();
std::vector<double> coordination4(nn);
std::vector<double> coordination2(nn);
std::vector<double> CVvalue_edge_temp(nn);
std::vector<double> CVvalue_face_temp(nn);
std::vector<double> CVvalue_corner_temp(nn);
std::vector<double> CVvalue_corner_angle(nn);
std::vector<double> Pb_ne(100);
int contact_mat[nn][12];
for (unsigned int k=0;k<nn;k+=1){
	for (unsigned int l=0;l<12;l+=1){
		contact_mat[k][l] = -1;
	}
}
int contact_mat_Pb[nn][12];
for (unsigned int k=0;k<nn;k+=1){
	for (unsigned int l=0;l<12;l+=1){
		contact_mat_Pb[k][l] = -1;
	}
}
 // Loop over center atoms
for(unsigned int i=0;i<nn;i+=1) { 	
        int k = 0;
        int p = 0;
	// Loop over around species
//	std::vector<unsigned> neighbors;
//	neighbors=nl->getNeighbors(i);
	for(unsigned int j=nn;j<(center_lista.size()+around_species1_lista.size());j+=1) {   
		double dfunc;
		Vector distance;
		Vector distance1;
		Vector distance2;
                double coordination1=0.;
		// i different from j
		//		if (atomsToRequest[i].serial()!=atomsToRequest[neighbors[j]].serial()) { 
		distance = pbcDistance(getPosition(i),getPosition(j));
		//		if ( (neighbors[j]>=center_lista.size()) && (neighbors[j]<(center_lista.size()+around_species1_lista.size())) ) {
		// Atom belongs to around species 1
		coordination1 = switchingFunction1.calculateSqr(distance.modulo2(),dfunc);
		coordination2[i] += switchingFunction1.calculateSqr(distance.modulo2(),dfunc);
	//	if (distance.modulo2() < 0.40){
		   if (coordination1 > 0.5){
			   contact_mat[i][k] = j;
			   k = k +1;
		}
	}
	if (k > 6){
		for(unsigned int l=0;l<12;l+=1) {
			if (contact_mat[i][l]!=-1) {
				for(unsigned int m=l;m<12-l;m+=1) {
					if (contact_mat[i][m+1]!=-1){
						Vector distance1 = pbcDistance(getPosition(i),getPosition(contact_mat[i][l]));
						Vector distance2 = pbcDistance(getPosition(i),getPosition(contact_mat[i][m+1])); 
						if(distance1.modulo2()>distance2.modulo2())
						{
							int temp=contact_mat[i][l];
							contact_mat[i][l]=contact_mat[i][m+1];
							contact_mat[i][m+1]=temp;
						}
					}
				}
			}
		}
	}
		for(unsigned int j=(center_lista.size()+around_species1_lista.size());j<(center_lista.size()+around_species1_lista.size()+around_species2_lista.size());j+=1) {   
		        double dfunc1;
			Vector distance1;
                        double coordination2=0.;
			distance1 = pbcDistance(getPosition(i),getPosition(j));
		        if ((distance1.modulo2()!= 0)) {
			coordination2 = switchingFunction2.calculateSqr(distance1.modulo2(),dfunc1);
			coordination4[i] += switchingFunction2.calculateSqr(distance1.modulo2(),dfunc1);
			   if (coordination2 > 0.5){
			   contact_mat_Pb[i][p] = j;
		           p = p +1;
	                } 
		     }
		}

	if (p > 6){
		for(unsigned int l=0;l<12;l+=1) {
			if (contact_mat_Pb[i][l]!=-1) {
				for(unsigned int m=l;m<12-l;m+=1) {
					if (contact_mat_Pb[i][m+1]!=-1){
						Vector distance1 = pbcDistance(getPosition(i),getPosition(contact_mat_Pb[i][l]));
						Vector distance2 = pbcDistance(getPosition(i),getPosition(contact_mat_Pb[i][m+1])); 
						if(distance1.modulo2()>distance2.modulo2())
						{
							int temp=contact_mat_Pb[i][l];
							contact_mat_Pb[i][l]=contact_mat_Pb[i][m+1];
							contact_mat_Pb[i][m+1]=temp;
						}
					}
				}
			}
		}
	}
             
}
//outputFile.printf("%f %f %f %f %f %f %f %f %f\n", getBox()[0][0], getBox()[0][1], getBox()[0][2], getBox()[1][0], getBox()[1][1], getBox()[1][2],getBox()[2][0], getBox()[2][1], getBox()[2][2]);
//outputFile.printField(nn);
//outputFile.printField("\n");
double corner = 0;
double face = 0;
double edge = 0;
double e1 = 0;
double e2 = 0;
double e3 = 0;
double e4 = 0;
double e5 = 0;
double angle = 0;
int ne = 0;
string some("Pb");
for(unsigned int i=0;i<nn;i+=1){
	if (coordination4[i] > 0.0){
		for (unsigned int k=0;k<nn;k+=1){
			if ((i!=k) && (coordination2[i] > 5.0) && (coordination2[k] > 5.0)){ 
				int temp = 0;
				int temp1 = 0;
				for(unsigned int j=0;j<6;j+=1){
					for(unsigned int m=0;m<6;m+=1){
						if ((contact_mat[i][j]!=-1) && (contact_mat[k][m]!=-1) && (contact_mat[i][j] == contact_mat[k][m]) ) {
							temp += 1;
						        double alpha;
						    	Vector distancei; 
						        Vector distancej;
						        distancei=pbcDistance(getPosition(i),getPosition(contact_mat[i][j]));
						        distancej=pbcDistance(getPosition(k),getPosition(contact_mat[i][j]));
							alpha = dotprod(distancei,distancej)/sqrt(norm2(distancei)*norm2(distancej));
                                                         if(alpha>=1.0){ 
                                                           alpha=0.99999;
                                                         }else if(alpha<=-1.0){
                                                           alpha=-0.99999;
                                                         }
                                                           double phi;
                                                           phi=acos(alpha);
							   if ((phi>=2.4)){
								   temp1 += 1;
							   }
						}
					}
				}
				if (temp1==1){
					CVvalue_corner_angle[i] += 1;
				}
				if (temp>=3){
					CVvalue_face_temp[i] += 1;
				} else if (temp == 2){
					CVvalue_edge_temp[i] += 1;
				} else if (temp == 1){
					CVvalue_corner_temp[i] += 1;  
				}
				//	log.printf("i, j %d %d\n",i, contact_mat[i][j]);
			}
		
		}

	}
}
for (int i=0;i<nn;i+=1){
		if (CVvalue_face_temp[i] >= 1){
		        face += 1;
		}
//	        if (CVvalue_edge_temp[i] >= 1){
//			e1 += 1;
//		}
//	        if (CVvalue_edge_temp[i] >= 2){
//			e2 += 1;
//		}
//	        if (CVvalue_edge_temp[i] >= 3){
//			e3 += 1;
//		}
//	        if (CVvalue_edge_temp[i] >= 4){
//			e4 += 1;
//		}
//	        if (CVvalue_edge_temp[i] >= 5){
//			e5 += 1;
//		}
		if (CVvalue_corner_temp[i] > 5.0){
		        corner += 1;
		}
	        if ((CVvalue_corner_angle[i] > 5.0)){
			 	angle += 1;
		                Pb_ne[ne]=i;
		                ne=ne+1;
		}
//outputFile.printf("%d Pb %f %f %f %f %f %f %f \n", i, getPosition(i)[0], getPosition(i)[1], getPosition(i)[2], CVvalue_face_temp[i], CVvalue_edge_temp[i], CVvalue_corner_temp[i], CVvalue_corner_angle[i]);
}

//outputFile.printf("Index, x, y, z, face, edge, corner %d %f %f %f %f %f %f \n", i, getPosition(i)[0], getPosition(i)[1], getPosition(i)[2], CVvalue_face, CVvalue_edge, CVvalue_corner);
//outputFile.printf("Pb %f %f %f %f %f %f %f \n", getPosition(i)[0], getPosition(i)[1], getPosition(i)[2], CVvalue_face_temp[i], CVvalue_edge_temp[i], CVvalue_corner_temp[i], CVvalue_corner_angle[i]);
//outputFile.printField(Pb).printField("x",getPosition(i)[0]).printField("y",getPosition(i)[1]).printField("z",getPosition(i)[2]).printField("face",CVvalue_face_temp[i]).printField("edge",CVvalue_edge_temp[i]).printField("corner",CVvalue_corner_temp[i]).printField();

//double surf = 0;
//for (int i=0;i<(ne-1);i+=1){
//       int m = Pb_ne[i];
//	for (int j=i+1;j<ne;j+=1){
//		int p = Pb_ne[j];
//		for (int k=0;k<6;k+=1){
//			for (int l=0;l<6;l+=1){
//                 		if ((contact_mat_Pb[m][k]!=-1) && (contact_mat_Pb[p][l]!=-1) && (contact_mat_Pb[p][l] == contact_mat_Pb[m][k])){
//		         	contact_mat_Pb[p][l] = -1;
//				}
//		}
//	}
//     }
//}

//for (int i=0;i<ne;i+=1) {
//      int m = Pb_ne[i];
////      outputFile.printf("Pb %f %f %f %d\n", getPosition(m)[0], getPosition(m)[1], getPosition(m)[2], 1);
//      int k = 0;      
//	while (k<6) {
////	        log.printf("i, j %d %d\n",m, (contact_mat_Pb[m][k]));
//		for (int j=0;j<ne;j+=1){
//			Vector distancek;
//			distancek=pbcDistance(getPosition(contact_mat_Pb[m][k]),getPosition(Pb_ne[j]));
////		        if (contact_mat_Pb[m][k] == (Pb_ne[j]+800)) {
//		        if (distancek.modulo2()==0) {
//				goto skip_loop;
//			}
//		}
//	   if (contact_mat_Pb[m][k]!=-1) {
//		surf += 1;
//		double neigh = contact_mat_Pb[m][k]; 
////                outputFile.printf("Pb %f %f %f %d\n", getPosition(neigh)[0], getPosition(neigh)[1], getPosition(neigh)[2], 2);
//	   }
//        skip_loop:
//	   k = k+1;
//	}
//}
//log.printf("face, e1, e2, e3, e4, e5, corner, angle, surface %f %f %f %f %f %f %f %f %f\n", face, e1, e2, e3, e4, e5, corner, angle, surf);
//outputFile.printf("\n");
//outputFile.close();
//if(!serial){
//   comm.Sum(CVvalue);
// }

setValue(corner);
}
double Oct::norm2(Vector vect){           /// CALCULATE THE NORM OF A VECTOR
  return (vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]);
}

double Oct::dotprod(Vector vect1,Vector vect2){           /// CALCULATE THE SCALAR PRODUCT OF TWO VECTORS
  return(vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2]);
}

}
}
