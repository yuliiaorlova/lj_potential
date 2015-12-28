#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>

using namespace std;

typedef struct
{
    float x, y, z;
} atom_position;

atom_position *positions = NULL;
//atom_position **F = NULL;
atom_position *F_res = NULL;
//float **En = NULL;
float *En_res=NULL;
int k=0, i=0,j=0;

float E=0, S=0, E4=0, E24=0, S2=0, S6=0, S12=0, S122=0;
string input_file_name, output_file_name;
int number_of_atoms;


void read_file(string input_file_name){
    
    ifstream input_file;
    string input[600];
  
    input_file.open(input_file_name.c_str());
    if(input_file.is_open()==0)
    {
        cout<<"error: could not create output file!"<<endl;
        exit(1);}
    
    while (!input_file.eof()) {
        getline(input_file, input[i]);
        i++;}
    k=i;
    number_of_atoms=atoi(input[0].c_str());
    cout<<"number of atoms="<< number_of_atoms<<endl;

    string *name_of_element;
    
    
    positions = new atom_position[number_of_atoms];
    F_res = new atom_position[number_of_atoms];
    En_res = new float[number_of_atoms];
    name_of_element = new string [number_of_atoms];
    
    std::istringstream iss;
    for (i=0; i<k-2; i++){
        iss.str(input[i+2]);
        iss >> name_of_element[i];
        iss >> positions[i].x;
        iss >> positions[i].y;
        iss >> positions[i].z;
        iss.clear();
    }

}


void calculate_forces(){
   float rij2=0, rij_2=0, rij_6=0, rij_12=0, diff_x=0, diff_y=0, diff_z=0, alpha=0;
    float Fx=0, Fy=0, Fz=0, Energ=0;
    for (i=0; i<number_of_atoms; i++)
    { for (j=0; j<number_of_atoms; j++)
    {
            if (i!=j)
            {
                
                diff_x=positions[i].x-positions[j].x;
                diff_y=positions[i].y-positions[j].y;
                diff_z=positions[i].z-positions[j].z;
                rij2=diff_x*diff_x+diff_y*diff_y+diff_z*diff_z;
                rij_2=1/rij2;
                rij_6=rij_2*rij_2*rij_2;
                rij_12=rij_6*rij_6;
                alpha=E24*rij_2*(S122*rij_12-S6*rij_6);
                
                Fx=alpha*diff_x;
                Fy=alpha*diff_y;
                Fz=alpha*diff_z;
                Energ=E4*(S12*rij_12-S6*rij_6);

            }
            else {
                
                Fx=0;
                Fy=0;
                Fz=0;
                Energ=0;
                
            }
            En_res[i]+=(0.5)*Energ;
            F_res[i].x+=Fx;
            F_res[i].y+=Fy;
            F_res[i].z+=Fz;
        }

    }
}


void write_file(string output_file_name){
    
    ofstream output_file;
    output_file.open(output_file_name.c_str());
    if(output_file.is_open()==0)
    {
        cout<<"error: could not create output file!"<<endl;
        exit(1);
    }
    for (i=0; i<number_of_atoms; i++){
        output_file<<F_res[i].x<<" "<<F_res[i].y<<" "<<F_res[i].z<<" "<<En_res[i]<<endl;
    }
    output_file.close();
}


int main(int argc, const char * argv[]) {
    clock_t t;
    t=clock();
    
    if (argc!=5){
        std::cout<<"error: wrong invocation!"<<endl;
        std::cout<<"try with:"<<endl;
        std::cout<<"potential <float E>, <float S>, <string input_file_name>, <string output_file_name>"<<endl;
        return 0;
    }
    
    E=atof(argv[1]);
    S=atof(argv[2]);
    input_file_name=argv[3];
    output_file_name=argv[4];
    
    cout<<"E="<<E<<endl;
    cout<<"S="<<S<<endl;
    cout<<"input_file_name="<<input_file_name<<endl;
    cout<<"output_file_name="<<output_file_name<<endl;
    
    //comoute constants
    E4=E*4;
    E24=E*24;
    S2=S*S;
    S6=S2*S2*S2;
    S12=S6*S6;
    S122=2*S12;
    
    
    //read_from_file
    read_file(input_file_name);

    //computation
     calculate_forces();
    
    //put data in file
    write_file(output_file_name);
    t=clock()-t;
    printf("it took me %lu clicks (%f seconds).\n", t, ((float)t)/CLOCKS_PER_SEC);
    //clean memory in the end of the program
    delete [] positions;
    delete [] F_res;
    delete [] En_res;
    
    
    return 0;
}
