#include<iostream>
#include<stack>
#include<vector>
#include<iomanip>
#include<fstream>
#include<ostream>
#include<math.h>
#include<ctime>
#include <cmath>
#include "bitmap_image.hpp"
using namespace std;
#define N 4
#define M 1
#define pi 2.0*acos(0.0)
#define INF numeric_limits<double>::infinity()

double eyeX, eyeY, eyeZ;
double lookX, lookY, lookZ;
double upX, upY, upZ;
double fovY, aspectRatio, near, far;
double** main_matrix;
double angle;
double ax,ay,az;
int counter = 0;
int Screen_Width,Screen_Height;
double left_limit_of_X,right_limit_of_X;
double bottom_limit_of_Y,top_limit_of_Y;
double front_limits_of_Z,rear_limits_of_Z;
double dx,dy;
double Top_Y,Left_X,Bottom_Y,Right_X;


struct point
{
	double x,y,z;

};
struct Color{

    int r,g,b;
};

struct Triangle
{
	struct point points[3];
    struct Color color;

};

//vector<vector<double> > two_d_vector;
stack<double**>s_tack;


void multiply_4x4(double** A,double** B,double** res)
{
	int i, j, k;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
            res[i][j] = 0;
			for (k = 0; k < N; k++)
				res[i][j] += A[i][k] * B[k][j];
		}
	}
}

void multiply_4x1(double** A,double** B,double** res)
{
	int i, j, k;
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
            res[i][j] = 0;
			for (k = 0; k < N; k++)
				res[i][j] += A[i][k] * B[k][j];
		}
	}
}

int compareStrings(string s1, string s2) {
   int compare = s1.compare(s2);
   if (compare != 0)
      return 0;
   else if(compare == 0)
      return 1;
}
void indentityMatrix(double** E){

    E[0][0] = 1;
    E[0][1] = 0;
    E[0][2] = 0;
    E[0][3] = 0;
    E[1][0] = 0;
    E[1][1] = 1;
    E[1][2] = 0;
    E[1][3] = 0;
    E[2][0] = 0;
    E[2][1] = 0;
    E[2][2] = 1;
    E[2][3] = 0;
    E[3][0] = 0;
    E[3][1] = 0;
    E[3][2] = 0;
    E[3][3] = 1;



}
void translateMatrix(double a,double b,double c,double **F){

    double** mat = new double*[N];
    for(int i = 0; i < N; ++i)
        mat[i] = new double[N];

    indentityMatrix(mat);
    mat[0][3]=a;
    mat[1][3]=b;
    mat[2][3]=c;
    int i,j;
    for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
            F[i][j]=mat[i][j];
		}
    }



}
void Normalize(double C[3],double D[3]){

        double temp;
        temp = sqrt(pow(C[0],2)+pow(C[1],2)+pow(C[2],2));
        D[0]=C[0]/temp;
        D[1]=C[1]/temp;
        D[2]=C[2]/temp;

}

double dot_product(double vector_a[], double vector_b[]) {
   double product = 0;
   for (int i = 0; i < 3; i++)
   product = product + vector_a[i] * vector_b[i];
   return product;
}

void cross_product(double vector_a[], double vector_b[], double temp[]) {
   temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
   temp[1] = -(vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0]);
   temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
   if(temp[0]==-0)
   {
       temp[0]=0;
   }
   if(temp[1]==-0)
   {
       temp[1]=0;
   }
   if(temp[2]==-0)
   {
       temp[2]=0;
   }

}

void viewTransformation(double **V){

    double look[3],eye[3],up[3],result[3],rW[3];
    double l[3],u[3],r[3];

    double** T = new double*[N];
    for(int i = 0; i < N; ++i)
        T[i] = new double[N];

    double** R = new double*[N];
    for(int i = 0; i < N; ++i)
        R[i] = new double[N];

    look[0]=lookX,look[1]=lookY,look[2]=lookZ;
    eye[0]=eyeX,eye[1]=eyeY,eye[2]=eyeZ;
    up[0]=upX,up[1]=upY,up[2]=upZ;

    for(int i=0;i<3;i++){
        l[i]=look[i]-eye[i];
    }

    Normalize(l,result);

    for(int i=0;i<3;i++){
        l[i]=result[i];
    }

    cross_product(l,up,r);

    Normalize(r,rW);

    for(int i=0;i<3;i++){
        r[i]=rW[i];
    }

    cross_product(r,l,u);

    translateMatrix(-eyeX,-eyeY,-eyeZ,T);

    R[0][0] = r[0];
    R[0][1] = r[1];
    R[0][2] = r[2];
    R[0][3] = 0;
    R[1][0] = u[0];
    R[1][1] = u[1];
    R[1][2] = u[2];
    R[1][3] = 0;
    R[2][0] = -l[0];
    R[2][1] = -l[1];
    R[2][2] = -l[2];
    R[2][3] = 0;
    R[3][0] = 0;
    R[3][1] = 0;
    R[3][2] = 0;
    R[3][3] = 1;

    double** V_ans = new double*[N];
    for(int i = 0; i < N; ++i)
        V_ans[i] = new double[N];

    multiply_4x4(R,T,V_ans);

    for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
                V[i][j]=V_ans[i][j];
		}
    }

}
void getC(double angle,double a[3],double x[3],double c[3])
{
    double arg = angle * (pi / 180.0);
    double s_ = sin( arg );
    double c_ = cos( arg );
    double temp[3];
    double temp_cross[3];
    double temp2[3],temp3[3];
    double t;

    for(int i=0;i<3;i++){
        temp[i]=c_*x[i];
    }

    double dot = dot_product(a,x);
    cross_product(a,x,temp_cross);

    t = (1-c_)*dot;

    for(int i=0;i<3;i++){
        temp2[i]=t*a[i];
    }

    for(int i=0;i<3;i++){
        temp3[i]=s_*temp_cross[i];
    }

    for(int i=0;i<3;i++){
        c[i]=temp[i]+temp2[i];
    }
    for(int i=0;i<3;i++){
        c[i]+=temp3[i];
    }

}

void projectionTransformation(double** proj_mat)
{
    double arg1 = (fovY/2) * (pi / 180.0);
    double t_1 = tan( arg1 );

    double fovX = fovY * aspectRatio;

    double arg2 = (fovX/2) * (pi / 180.0);
    double t_2 = tan( arg2 );

    double t = near * t_1;
    double r = near * t_2;
    double temp =-((far+near)/(far-near));
    double temp1 =-(((2*far)*near)/(far-near));


    proj_mat[0][0] = near/r;
    proj_mat[0][1] = 0;
    proj_mat[0][2] = 0;
    proj_mat[0][3] = 0;
    proj_mat[1][0] = 0;
    proj_mat[1][1] = near/t;
    proj_mat[1][2] = 0;
    proj_mat[1][3] = 0;
    proj_mat[2][0] = 0;
    proj_mat[2][1] = 0;
    proj_mat[2][2] = temp;
    proj_mat[2][3] = temp1;
    proj_mat[3][0] = 0;
    proj_mat[3][1] = 0;
    proj_mat[3][2] = -1;
    proj_mat[3][3] = 0;


}
double interpolation(double ys,double y1,double y3,double z1,double z3){

  double answer;
  if((y3-y1)== 0.0){
        answer = z1;
  }
  else if((ys-y1)== 0.0 or (z3-z1)== 0.0){
        answer = z1;
  }
  else{
    double ans1,ans2,ans3,ans4;

     ans1 = z3-z1;
     ans2 = y3-y1;
     ans3 = ys-y1;
     ans4 = ans1*ans3;
     answer = ans4/ans2;
     answer = answer+z1;

  }

    return answer;


}
bool inRange(double l, double h, double x)
{
double low,high;
    if(l<=h){
        low = l,high =h;
        return (low <= x && x <= high);
    }
    else{
        low = h,high =l;
        return (low <= x && x <= high);
    }

}

int main()
{
    //double res[N][N];
    main_matrix = new double*[N];
    for(int i = 0; i < N; ++i)
        main_matrix[i] = new double[N];

    main_matrix[0][0] = 1;
    main_matrix[0][1] = 0;
    main_matrix[0][2] = 0;
    main_matrix[0][3] = 0;
    main_matrix[1][0] = 0;
    main_matrix[1][1] = 1;
    main_matrix[1][2] = 0;
    main_matrix[1][3] = 0;
    main_matrix[2][0] = 0;
    main_matrix[2][1] = 0;
    main_matrix[2][2] = 1;
    main_matrix[2][3] = 0;
    main_matrix[3][0] = 0;
    main_matrix[3][1] = 0;
    main_matrix[3][2] = 0;
    main_matrix[3][3] = 1;


    /*res[N][N] ={ { 0, 0, 0, 0 },
               { 0, 0, 0, 0 },
               { 0, 0, 0, 0 },
               { 0, 0, 0, 0 } };*/

    ifstream fin;
    fin.open("scene.txt",ios::in);


    fstream stage1,stage2,stage3;
    stage1.open("stage1.txt",std::ios_base::app);
    stage2.open("stage2.txt",std::ios_base::app);
    stage3.open("stage3.txt",std::ios_base::app);

    fin>>eyeX>>eyeY>>eyeZ;
    fin>>lookX>>lookY>>lookZ;
    fin>>upX>>upY>>upZ;
    fin>>fovY>>aspectRatio>>near>>far;

    double** matri_x = new double*[N];
    for(int i = 0; i < N; ++i)
        matri_x[i] = new double[N];

    viewTransformation(matri_x);

    double** projectionMatrix = new double*[N];
    for(int i = 0; i < N; ++i)
        projectionMatrix[i] = new double[N];

    projectionTransformation(projectionMatrix);

    while(1)
    {
        string choice;
        fin>>choice;

        if(compareStrings(choice,"push")==1){

            double** tempo = new double*[N];
            for(int i = 0; i < N; ++i)
               tempo[i] = new double[N];

            //cout<<"push"<<endl;
            for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        tempo[i][j]=main_matrix[i][j];
                        //cout<<tempo[i][j]<<" ";
                    }
                    //cout<<" "<<endl;
              }

            s_tack.push(tempo);

        }
        else if(compareStrings(choice,"pop")==1){

              /*double** tempo = new double*[N];
                for(int i = 0; i < N; ++i)
               tempo[i] = new double[N];*/
                /*cout<<"before popped"<<endl;
               for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        cout<<main_matrix[i][j]<<" ";
                    }cout<<" "<<endl;
              }*/

               main_matrix =s_tack.top();
                s_tack.pop();

               /* cout<<"after popped"<<endl;
               for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        cout<<main_matrix[i][j]<<" ";
                    }cout<<" "<<endl;
              }*/




        }

        else if(compareStrings(choice,"scale")==1){

            double sx,sy,sz;
            fin>>sx>>sy>>sz;

            double** temp = new double*[N];
            for(int i = 0; i < N; ++i)
               temp[i] = new double[N];

            double** scale = new double*[N];
            for(int i = 0; i < N; ++i)
               scale[i] = new double[N];

            temp[0][0]=sx;
            temp[0][1]=0;
            temp[0][2]=0;
            temp[0][3]=0;
            temp[1][0]=0;
            temp[1][1]=sy;
            temp[1][2]=0;
            temp[1][3]=0;
            temp[2][0]=0;
            temp[2][1]=0;
            temp[2][2]=sz;
            temp[2][3]=0;
            temp[3][0]=0;
            temp[3][1]=0;
            temp[3][2]=0;
            temp[3][3]=1;

            multiply_4x4(main_matrix,temp,scale);

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                        main_matrix[i][j]=scale[i][j];
                }
            }

        }
        else if(choice=="translate"){

            double tx,ty,tz;
            fin>>tx>>ty>>tz;

            double** temp = new double*[N];
            for(int i = 0; i < N; ++i)
               temp[i] = new double[N];

            double** trans = new double*[N];
            for(int i = 0; i < N; ++i)
               trans[i] = new double[N];


            temp[0][0]=1;
            temp[0][1]=0;
            temp[0][2]=0;
            temp[0][3]=tx;
            temp[1][0]=0;
            temp[1][1]=1;
            temp[1][2]=0;
            temp[1][3]=ty;
            temp[2][0]=0;
            temp[2][1]=0;
            temp[2][2]=1;
            temp[2][3]=tz;
            temp[3][0]=0;
            temp[3][1]=0;
            temp[3][2]=0;
            temp[3][3]=1;

            multiply_4x4(main_matrix,temp,trans);

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                        main_matrix[i][j]=trans[i][j];
                        //printf("%1f \n",main_matrix[i][j]);
                }
            }


        }
        else if(compareStrings(choice,"rotate")==1){

            fin>>angle>>ax>>ay>>az;
            double i_[3],j_[3],k[3];
            double c1[3],c2[3],c3[3],a_main[3],a[3];

            a_main[0]=ax;
            a_main[1]=ay;
            a_main[2]=az;

            Normalize(a_main,a);

            i_[0]=1;
            i_[1]=0;
            i_[2]=0;

            j_[0]=0;
            j_[1]=1;
            j_[2]=0;

            k[0]=0;
            k[1]=0;
            k[2]=1;

            getC(angle,a,i_,c1);
            getC(angle,a,j_,c2);
            getC(angle,a,k,c3);

            double** rotationMatrix= new double*[N];
            for(int i = 0; i < N; ++i)
               rotationMatrix[i] = new double[N];

            rotationMatrix[0][0]=c1[0];
            rotationMatrix[0][1]=c2[0];
            rotationMatrix[0][2]=c3[0];
            rotationMatrix[0][3]=0;
            rotationMatrix[1][0]=c1[1];
            rotationMatrix[1][1]=c2[1];
            rotationMatrix[1][2]=c3[1];
            rotationMatrix[1][3]=0;
            rotationMatrix[2][0]=c1[2];
            rotationMatrix[2][1]=c2[2];
            rotationMatrix[2][2]=c3[2];
            rotationMatrix[2][3]=0;
            rotationMatrix[3][0]=0;
            rotationMatrix[3][1]=0;
            rotationMatrix[3][2]=0;
            rotationMatrix[3][3]=1;

            double** ans = new double*[N];
            for(int i = 0; i < N; ++i)
               ans[i] = new double[N];

            multiply_4x4(main_matrix,rotationMatrix,ans);

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    main_matrix[i][j]=ans[i][j];
                }
            }


        }
        else if(compareStrings(choice,"triangle")==1){

            counter +=1;
            //cout<<"counter"<<endl;

            double** tri_angle1 = new double*[N];
            for(int i = 0; i < N; ++i)
               tri_angle1[i] = new double[N];

            double** tri_angle2 = new double*[N];
            for(int i = 0; i < N; ++i)
               tri_angle2[i] = new double[N];

            double** tri_angle3 = new double*[N];
            for(int i = 0; i < N; ++i)
               tri_angle3[i] = new double[N];


            fin>>tri_angle1[0][0]>>tri_angle1[1][0]>>tri_angle1[2][0];
            fin>>tri_angle2[0][0]>>tri_angle2[1][0]>>tri_angle2[2][0];
            fin>>tri_angle3[0][0]>>tri_angle3[1][0]>>tri_angle3[2][0];


            tri_angle1[3][0]=1;
            tri_angle2[3][0]=1;
            tri_angle3[3][0]=1;

            double** t1 = new double*[N];
            for(int i = 0; i < N; ++i)
               t1[i] = new double[N];

            double** t2 = new double*[N];
            for(int i = 0; i < N; ++i)
               t2[i] = new double[N];

            double** t3 = new double*[N];
            for(int i = 0; i < N; ++i)
               t3[i] = new double[N];


            int i1,j1,i2,j2,i3,j3;

            multiply_4x1(main_matrix,tri_angle1,t1);
            multiply_4x1(main_matrix,tri_angle2,t2);
            multiply_4x1(main_matrix,tri_angle3,t3);



            if(t1[3][0]!=1.000000){
                t1[0][0]=t1[0][0]/t1[3][0];
                t1[1][0]=t1[1][0]/t1[3][0];
                t1[2][0]=t1[2][0]/t1[3][0];
                //t1[0][3]=t1[0][3]/t1[0][3];
            }

            if(t2[3][0]!=1.000000){
                t2[0][0]=t2[0][0]/t2[3][0];
                t2[1][0]=t2[1][0]/t2[3][0];
                t2[2][0]=t2[2][0]/t2[3][0];
                //t1[0][3]=t1[0][3]/t1[0][3];
            }

            if(t3[3][0]!=1.000000){
                t3[0][0]=t3[0][0]/t3[3][0];
                t3[1][0]=t3[1][0]/t3[3][0];
                t3[2][0]=t3[2][0]/t3[3][0];
                //t1[0][3]=t1[0][3]/t1[0][3];
            }


            for (i1 = 0; i1 < N-1; i1++) {
            for (j1 = 0; j1 < M; j1++)
			stage1<<std::fixed<<setprecision(7)<<t1[i1][j1]<<" ";
            }
            stage1<<" "<<endl;

            for (i2 = 0; i2 < N-1; i2++) {
            for (j2 = 0; j2 < M; j2++)
			stage1<<std::fixed<<setprecision(7)<<t2[i2][j2]<<" ";
            }
            stage1<<" "<<endl;

            for (i3 = 0; i3 < N-1; i3++) {
            for (j3 = 0; j3 < M; j3++)
			stage1<<std::fixed<<setprecision(7)<<t3[i3][j3]<<" ";
            }
            stage1<<" "<<endl;
            stage1<<" "<<endl;


            double** p1 = new double*[N];
            for(int i = 0; i < N; ++i)
               p1[i] = new double[N];

            double** p2 = new double*[N];
            for(int i = 0; i < N; ++i)
               p2[i] = new double[N];

            double** p3 = new double*[N];
            for(int i = 0; i < N; ++i)
               p3[i] = new double[N];

            multiply_4x1(matri_x,t1,p1);
            multiply_4x1(matri_x,t2,p2);
            multiply_4x1(matri_x,t3,p3);

            if(p1[3][0]!=1.000000){
                p1[0][0]=p1[0][0]/p1[3][0];
                p1[1][0]=p1[1][0]/p1[3][0];
                p1[2][0]=p1[2][0]/p1[3][0];
                //t1[0][3]=t1[0][3]/t1[0][3];
            }

            if(p2[3][0]!=1.000000){
                p2[0][0]=p2[0][0]/p2[3][0];
                p2[1][0]=p2[1][0]/p2[3][0];
                p2[2][0]=p2[2][0]/p2[3][0];
                //t1[0][3]=t1[0][3]/t1[0][3];
            }

            if(p3[3][0]!=1.000000){
                p3[0][0]=p3[0][0]/p3[3][0];
                p3[1][0]=p3[1][0]/p3[3][0];
                p3[2][0]=p3[2][0]/p3[3][0];
                //t1[0][3]=t1[0][3]/t1[0][3];
            }

            for (int i = 0; i< N-1; i++) {
            for (int j = 0; j < M; j++)
			stage2<<std::fixed<<setprecision(7)<<p1[i][j]<<" ";
            }
            stage2<<" "<<endl;

            for (int i = 0; i < N-1; i++) {
            for (int j = 0; j < M; j++)
			stage2<<std::fixed<<setprecision(7)<<p2[i][j]<<" ";
            }
            stage2<<" "<<endl;

            for (int i = 0; i < N-1; i++) {
            for (int j = 0; j < M; j++)
			stage2<<std::fixed<<setprecision(7)<<p3[i][j]<<" ";
            }
            stage2<<" "<<endl;
            stage2<<" "<<endl;


            //eikhan theke projection

            double** n1 = new double*[N];
            for(int i = 0; i < N; ++i)
               n1[i] = new double[N];

            double** n2 = new double*[N];
            for(int i = 0; i < N; ++i)
               n2[i] = new double[N];

            double** n3 = new double*[N];
            for(int i = 0; i < N; ++i)
               n3[i] = new double[N];

            multiply_4x1(projectionMatrix,p1,n1);
            multiply_4x1(projectionMatrix,p2,n2);
            multiply_4x1(projectionMatrix,p3,n3);

            if(n1[3][0]!=1.000000){
                n1[0][0]=n1[0][0]/n1[3][0];
                n1[1][0]=n1[1][0]/n1[3][0];
                n1[2][0]=n1[2][0]/n1[3][0];
                //t1[0][3]=t1[0][3]/t1[0][3];
            }

            if(n2[3][0]!=1.000000){
                n2[0][0]=n2[0][0]/n2[3][0];
                n2[1][0]=n2[1][0]/n2[3][0];
                n2[2][0]=n2[2][0]/n2[3][0];
                //t1[0][3]=t1[0][3]/t1[0][3];
            }

            if(n3[3][0]!=1.000000){
                n3[0][0]=n3[0][0]/n3[3][0];
                n3[1][0]=n3[1][0]/n3[3][0];
                n3[2][0]=n3[2][0]/n3[3][0];
                //t1[0][3]=t1[0][3]/t1[0][3];
            }

            for (int i = 0; i< N-1; i++) {
            for (int j = 0; j < M; j++)
			stage3<<std::fixed<<setprecision(7)<<n1[i][j]<<" ";
            }
            stage3<<" "<<endl;

            for (int i = 0; i < N-1; i++) {
            for (int j = 0; j < M; j++)
			stage3<<std::fixed<<setprecision(7)<<n2[i][j]<<" ";
            }
            stage3<<" "<<endl;

            for (int i = 0; i < N-1; i++) {
            for (int j = 0; j < M; j++)
			stage3<<std::fixed<<setprecision(7)<<n3[i][j]<<" ";
            }
            stage3<<" "<<endl;
            stage3<<" "<<endl;


        }
        else if(compareStrings(choice,"end")==1){
                stage1.close();
                stage2.close();
                stage3.close();
                break;

        }
    }


        ifstream fstage3;
        fstage3.open("stage3.txt",ios::in);
        fstream z_buffer;
        z_buffer.open("z-buffer.txt",std::ios_base::app);
        ifstream config;
        config.open("config.txt",ios::in);
        config>>Screen_Width>>Screen_Height;
        config>>left_limit_of_X;
        config>>bottom_limit_of_Y;
        config>>front_limits_of_Z>>rear_limits_of_Z;

        bitmap_image image(Screen_Width,Screen_Height);
        //image.save_image("output.bmp");

        double** Z_buffer = new double*[Screen_Height];
        for(int i = 0; i < Screen_Height; ++i)
            Z_buffer[i] = new double[Screen_Width];


        double z_max = rear_limits_of_Z;

        for(int i=0;i<Screen_Height;i++){
            for(int j=0;j<Screen_Width;j++){
                Z_buffer[i][j] = z_max;
            }
        }

       Color** frameBuffer = new Color*[Screen_Height];
        for(int i=0; i<Screen_Height; i++) {
            frameBuffer[i] = new Color[Screen_Width];
        }

        for(int i=0; i<Screen_Height; i++) {
            for(int j=0; j<Screen_Width; j++) {
            frameBuffer[i][j].r = 0;
            frameBuffer[i][j].g = 0;
            frameBuffer[i][j].b = 0;
            }
        }

        right_limit_of_X=-(left_limit_of_X);
        top_limit_of_Y=-(bottom_limit_of_Y);

        dx = (right_limit_of_X - left_limit_of_X) / Screen_Width;
        dy = (top_limit_of_Y - bottom_limit_of_Y) / Screen_Height;

        Top_Y = top_limit_of_Y - (dy/2);
        Bottom_Y = bottom_limit_of_Y + (dy/2);
        Left_X = left_limit_of_X + (dx/2);
        Right_X = right_limit_of_X - (dx/2);

        int Top_scanline,Bottom_scanline;

        struct Triangle triangles_[counter];

        srand(time(0));

        for(int i=0;i<counter;i++){
           // cout<<"TRiangle start"<<endl;

            fstage3>>triangles_[i].points[0].x>>triangles_[i].points[0].y>>triangles_[i].points[0].z;
            fstage3>>triangles_[i].points[1].x>>triangles_[i].points[1].y>>triangles_[i].points[1].z;
            fstage3>>triangles_[i].points[2].x>>triangles_[i].points[2].y>>triangles_[i].points[2].z;




            triangles_[i].color.r= rand()%256;
            triangles_[i].color.g= rand()%256;
            triangles_[i].color.b= rand()%256;


            double max_y= max(triangles_[i].points[0].y,triangles_[i].points[1].y);
            max_y = max(max_y,triangles_[i].points[2].y);
            double min_y= min(triangles_[i].points[0].y,triangles_[i].points[1].y);
            min_y = min(min_y,triangles_[i].points[2].y);


            if(max_y >= Top_Y) {
            Top_scanline = 0;
            }
            else {
            Top_scanline = (int) round((Top_Y - max_y)/dy);
            }
            if(min_y <= Bottom_Y) {
                Bottom_scanline = Screen_Height - 1;
            }
            else {
            Bottom_scanline = Screen_Height - (1 + ((int) round((min_y - Bottom_Y)/dy)));
            }

            //cout<<"TOPSCAN :"<<Top_scanline<<"BOTTOMSCAN :"<<Bottom_scanline<<endl;

            for(int k = Top_scanline; k<= Bottom_scanline; k++){

                point p,q;
                double y_s = Top_Y-(k*dy);
                //cout<<"y_s"<<y_s<<endl;
                p.y = q.y = y_s;

                int tau1=0;
                int tau2=1;
                int tau3=0;
                int tau4=2;

                if(inRange(triangles_[i].points[tau1].y,triangles_[i].points[tau2].y,y_s)==1)
                {
                    if(inRange(triangles_[i].points[tau3].y,triangles_[i].points[tau4].y,y_s)==1){

                        continue;
                    }
                    else{
                        tau3 = 1;
                        tau4 = 2;
                    }
                }
                else{

                    tau1 = 1;
                    tau2 = 2;

                }
                /*cout<<triangles_[i].points[0].x<<" "<<triangles_[i].points[0].y<<" "<<triangles_[i].points[0].z<<endl;
                cout<<triangles_[i].points[1].x<<" "<<triangles_[i].points[1].y<<" "<<triangles_[i].points[1].z<<endl;
                cout<<triangles_[i].points[2].x<<" "<<triangles_[i].points[2].y<<" "<<triangles_[i].points[2].z<<endl;
                cout<<" "<<endl;

                cout << "taus " << tau1 << ' ' << tau2 << ' ' << tau3 << ' ' << tau4 << endl;*/
                p.x = interpolation(y_s,triangles_[i].points[tau1].y,triangles_[i].points[tau2].y,triangles_[i].points[tau1].x,triangles_[i].points[tau2].x);
                p.z = interpolation(y_s,triangles_[i].points[tau1].y,triangles_[i].points[tau2].y,triangles_[i].points[tau1].z,triangles_[i].points[tau2].z);

                q.x =interpolation(y_s,triangles_[i].points[tau3].y,triangles_[i].points[tau4].y,triangles_[i].points[tau3].x,triangles_[i].points[tau4].x);
                q.z = interpolation(y_s,triangles_[i].points[tau3].y,triangles_[i].points[tau4].y,triangles_[i].points[tau3].z,triangles_[i].points[tau4].z);

                //cout<<"p.x "<<p.x<<"p.z "<<p.z<<"q.x "<<q.x<<"q.z "<<q.z<<endl;

                double p_mid = (p.x - left_limit_of_X) / dx;
                double q_mid = (q.x - left_limit_of_X) / dx;

                int p_ = (int) round(p_mid);
                int q_ = (int) round(q_mid);
                int temp7= max(p_,q_);
                int temp8= min(p_,q_);
                p_ =temp8;
                q_ =temp7;
                //cout<<"p_"<<p_<<"q_"<<q_<<endl;
                //eikhane somossha.p_x ,q_x er value protibar 0 ashe.tai r dhukteche na loop a

                for(int j =p_;j<=q_;j++)//ei line a jhamela
                {
                    //cout<<"in the inner"<<endl;
                    double value_ = (j*dx)+Left_X;
                    double z_val = interpolation(value_,p.x,q.x,p.z,q.z);
                    if(z_val < front_limits_of_Z){
                        continue;
                    }
                    if(z_val > rear_limits_of_Z){
                        continue;
                    }
                    if(z_val<Z_buffer[k][j]){
                        Z_buffer[k][j]=z_val;

                        frameBuffer[k][j].r = triangles_[i].color.r;
                        frameBuffer[k][j].g = triangles_[i].color.g;
                        frameBuffer[k][j].b = triangles_[i].color.b;
                    }
                }


            }


            cout<<triangles_[i].points[0].x<<" "<<triangles_[i].points[0].y<<" "<<triangles_[i].points[0].z<<endl;
            cout<<triangles_[i].points[1].x<<" "<<triangles_[i].points[1].y<<" "<<triangles_[i].points[1].z<<endl;
            cout<<triangles_[i].points[2].x<<" "<<triangles_[i].points[2].y<<" "<<triangles_[i].points[2].z<<endl;
            cout<<" "<<endl;

            //cout<<"Triangle end"<<endl;
        }

        for (int l = 0; l< Screen_Height; l++) {
             for (int p = 0; p < Screen_Width; p++){
                    z_buffer<<std::fixed<<setprecision(7)<<Z_buffer[l][p]<<" ";
                    image.set_pixel(p,l, frameBuffer[l][p].r,frameBuffer[l][p].g,frameBuffer[l][p].b);
                 }
             }
        z_buffer<<" "<<endl;

        image.save_image("output.bmp");

       z_buffer.close();

       stage3.close();


    for(int i=0; i<N; i++) {
        delete[] main_matrix[i];
    }
    delete[] main_matrix;

    for(int i=0; i<Screen_Height; i++) {
        delete[] Z_buffer[i];
        delete[] frameBuffer[i];
    }
    delete[] Z_buffer;
    delete[] frameBuffer;





	return 0;
}


