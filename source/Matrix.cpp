#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

namespace mc1{
          //! Miscellanea container template class
          template <class type>
          class utilities{
                public:
                       //! e^'x' approximation with the first 'terms' terms of the series
                       type expMacLaurin(type x, unsigned long int terms){
                            type result=1;
                            type term=1;
                            for(int counter=1; counter<terms; counter++){
                                    term*=x;
                                    term/=counter;
                                    result+=term;
                                    }
                            return result;
                            }
                            
                       //! Absolute difference
                       type absDiff(type testval, type x){
                            return testval-x;
                            }
                            
                       //! Relative difference
                       type relDiff(type testval, type x){
                            return absDiff(testval,x)/testval;
                            }
                            
                       //! Machine epsilon
                       type epsilon(){
                            unsigned int d=0;
                            type eps=1;
                            while(1+eps>1){
                                           eps/=2;
                                           d+=1;
                                           }
                                           return d-1;
                            }
                        
                };
          
          //! Vector mask template class
          //! Enables the vector operations over arrays of data
          //! Warning: if the data is not embedded, the operations will be carried out over the original
          template <class type>
          class vector{
                private:
                        type *array; //! Pointer to data location
                        bool embedded; //! If true, 'array' is deallocated as the object is destroyed
                public:
                       int length; //! Length of the vector
                       
                       //! Builds a vector over the array 'a' of size 'l' [3]
                       vector(type *a, int l){
                                   array=a;
                                   length=l;
                                   embedded=false;
                                   }
                                   
                       //! Builds a vector of size 'l' with embedded uninitialized data [length+2]
                       vector(int l){
                                   array=(type *)calloc(l,sizeof(type));
                                   length=l;
                                   embedded=true;
                                  }
                                  
                       //! Destructor
                       ~vector(){
                                 if(embedded){
                                              free(array);
                                              }
                                 }
                                 
                       //! Prints to screen line the content of the vector [2*length+1]
                       void print(){
                            for(int i=0; i<length; i++){
                                    cout<<array[i]<<" ";
                                    }
                            cout<<endl;
                            }
                            
                       //! Makes the data embedded [1]
                       void embed(){
                            embedded=true;
                            }
                            
                       //! Makes the data not embedded [1]
                       void unembed(){
                            embedded=false;
                            }
                            
                       //! Checks the embedding of data [1]
                       bool isEmbedded(){
                            return embedded;
                            }
                            
                       //! Access operator [1]
                       type& operator[](int index){
                             return array[index];
                             }
                             
                       //! Copies the data into a pre-allocated array [length]
                       void copyIntoArray(type *target){
                            for(int i=0; i<length; i++){
                                    target[i]=array[i];
                                    }
                            }
                            
                       //! Returns a pointer to the address of a vector, equal to this, with embedded copy of data [2*length+4]
                       vector *clone(){
                              type *nn=(type *) calloc(length,sizeof(type));
                              copyIntoArray(nn);
                              vector *vv=new vector(nn,length);
                              vv->embed();
                              return vv;
                              }
                        
                        //! returns a stack copy of this vector      
                        vector stackClone(){
                            type *nn=(type *) calloc(length,sizeof(type));
                              copyIntoArray(nn);
                              vector vv(nn,length);
                              vv.embed();
                              return vv;
                        }
                              
                       //! Vector addition [2*length+2]
                       void add(vector *v){
                            int d=min(length,v->length);
                            for(int i=0; i<d; i++){
                                    array[i]+=(*v)[i];
                                    }
                            }
                            
                       //! Vector subtraction [2*length+2]
                       void sub(vector *v){
                            int d=min(length,v->length);
                            for(int i=0; i<d; i++){
                                    array[i]-=(*v)[i];
                                    }
                            }
                            
                       //! Vector inversion [2*length+1]
                       void inverse(){
                            int d=length;
                            for(int i=0; i<d; i++){
                                    array[i]-=array[i];
                                    }
                            }
                            
                       //! Standard external product [3*length+3]
                       void prod(type val){
                            for(int i=0; i<length; i++){
                                    array[i]*=val;
                                    }
                            }
                               
                       //! Returns a pointer to the address of a vector containing absolute differences with vector v
                       vector *absDiffs(vector *v){
                              type *nn=(type *) calloc(length,sizeof(type));
                              copyIntoArray(nn);
                              vector *vv=new vector(nn,length);
                              vv->embed();
                              vv->sub(v);
                              return vv;
                              }
                              
                       //! Returns a pointer to the address of a vector containing relative differences with vector v
                       vector *relDiffs(vector *v){
                              type *nn=(type *) calloc(length,sizeof(type));
                              copyIntoArray(nn);
                              vector *vv=new vector(nn,length);
                              vv->embed();
                              vv->sub(v);
                              for(int d=0; d<vv->length; d++){
                                      (*vv)[d]/=(*v)[d];
                                      }
                              return vv;
                              }
                                   
                       //! Standard Euclidean scalar product [3*length+3]
                       type prod(vector *v){
                            int d=min(length,v->length);
                            type result=0;
                            for(int i=0; i<d; i++){
                                    result+=array[i]*(*v)[i];
                                    }
                            return result;
                            }
                            
                       //! Resizes the vector, with zero padding in case of extension [3*length+6]
                       void resize(int newsize){
                            type *nn=calloc(newsize,sizeof(type));
                            for(int i=0; i<newsize; i++){
                                    nn[i]=0;
                                    }
                            length=min(length,newsize);
                            copyIntoArray(nn);
                            if(embedded) free(array);
                            length=newsize;
                            array=nn;
                            embed();
                            }
                            
                       //! Summation of all elements [length+1]
                       type summation(){
                            type res=0;
                            for(int i=0; i<length; i++){
                                    res+=array[i];
                                    }
                            return res;
                            }
                        
                       //! Vector sum
                       vector operator+(vector v){
                            vector vv=stackClone();
                            vv.add(&v);
                            return vv;
                        }
                        
                        //! external product
                        vector operator*(type val){
                            vector vv=stackClone();
                            vv.prod(val);
                            return vv;
                        }
                        
                        //! Vector difference
                       vector operator-(vector v){
                            vector vv=v*-1;
                            vv.add(this);
                            return vv;
                        }
                        
                        //! external division
                        vector operator/(type val){
                            vector vv=stackClone();
                            vv.prod(1/val);
                            return vv;
                        }
                        
                };
                
          //! 2d line equation container
          template <class type>
          class linearEquation2d{
                private:
                    type coefficients[3];
                    type sign;
                public:
                    linearEquation2d(type coeffX, type coeffY, type Kcoeff){
                        sign=0;
                        coefficients[0]=coeffX;
                        coefficients[1]=coeffY;
                        coefficients[2]=Kcoeff;
                        sign=+1;
                    }
                    linearEquation2d(mc1::vector<type> first, mc1::vector<type> second){
                        sign=0;
                        type x1=first[0];
                        type x2=second[0];
                        type y1=first[1];
                        type y2=second[1];
                        if(y1==y2){
                            type coeffX=0;
                            type coeffY=1;
                            type Kcoeff=-y1;
                            coefficients[0]=coeffX;
                            coefficients[1]=coeffY;
                            coefficients[2]=Kcoeff;
                            if(x1==x2) return;
                            if(x1<x2) sign=+1;
                            else sign=-1;
                        }
                        else{
                            type coeffX=1;
                            type coeffY=(x2-x1)/(y1-y2);
                            type Kcoeff=-(coeffY*y1+x1);
                            coefficients[0]=coeffX;
                            coefficients[1]=coeffY;
                            coefficients[2]=Kcoeff;
                            if(y1<y2) sign=+1;
                            else sign=-1;
                        }
                    }
                    type eval(type x, type y){
                        return sign*(coefficients[0]*x+coefficients[1]*y+coefficients[2]);
                    }
                    type eval(mc1::vector<type> v){
                        type x=v[0]; 
                        type y=v[1];
                        return sign*(coefficients[0]*x+coefficients[1]*y+coefficients[2]);
                    }
            };
                
          //! Matrix mask template class
          //! Enables the matrix operations over bidimensional arrays (type**) of data
          //! Warning: if the data is not embedded, the operations will be carried out over the original    
          template <class type>
          class matrix{
                private:
                        type **table; //! Pointer to data location
                        mutable bool embedded; //! If true, 'table' is deallocated as the object is destroyed
                public:
                        int rows; //! Rows of the matrix
                        int cols; //! Columns of the matrix
                        
                       //! Builds a matrix over a table 't' with 'r' rows and 'c' columns [4]
                       matrix(type **t, int r, int c){
                                   table=t;
                                   rows=r;
                                   cols=c;
                                   embedded=false;
                                   }
                       
                        //! Copy constructor
                       matrix(const matrix& matr){
                                    rows=matr.rows;
                                    cols=matr.cols;
                                   table=matr.table;
                                   matr.embedded=false;
                                   embedded=true;
                                   }
                                               
                       //! Builds an Hilbert matrrix of 'rr' rows and 'cc' columns[7*rows*cols+2]
                       matrix(int rr, int cc){
                                  table =(type **) calloc(rr,sizeof(type*));
                                  rows=rr;
                                  cols=cc;
                                  for(int r=0; r<rows; r++){
                                          table[r]=(type *) calloc(cols, sizeof(type));
                                          for(int c=0; c<cols; c++){
                                                       table[r][c]=((type) 1.0)/((type) r+c+1.0);
                                                  }
                                          }
                                  }
                       //! Builds a Givens rotation matrix of 'rr' rows and 'cc' columns[7*rows*cols+2]
                       matrix(int rc, int pos1, int pos2, type cos, type sin){
                                  table =(type **) calloc(rc,sizeof(type*));
                                  rows=rc;
                                  cols=rc;
                                  for(int r=0; r<rows; r++){
                                          table[r]=(type *) calloc(cols, sizeof(type));
                                          for(int c=0; c<cols; c++){
                                                       table[r][c]=(r==c?1:0);
                                                  }
                                          }
                                  if(pos1<0 || pos1>=rc) return;
                                  if(pos2<0 || pos2>=rc) return;
                                  table[pos1][pos1]=table[pos2][pos2]=cos;
                                  table[pos1][pos2]=-(table[pos2][pos1]=sin);
                                  }
                       //! Builds an Identity square matrix of 'rc' rows and 'rc' columns[7*rows*cols+2]
                       matrix(int rc){
                                  table =(type **) calloc(rc,sizeof(type*));
                                  rows=rc;
                                  cols=rc;
                                  for(int r=0; r<rows; r++){
                                          table[r]=(type *) calloc(cols, sizeof(type));
                                          for(int c=0; c<cols; c++){
                                                       table[r][c]=(r==c?1:0);
                                                  }
                                          }
                                  }          
                       //! Destructor
                       ~matrix(){
                                 if(embedded){
                                              for(int r=0; r<rows; r++){
                                                      free(table[r]);
                                                      }
                                              free(table);
                                              }
                                 }
                       
                       //! Embedding setters/getter
                       void embed(){
                            embedded=true;
                            }
                       void unembed(){
                            embedded=false;
                            }
                       bool isEmbedded(){
                            return embedded;
                            }
                       
                       //! Access operators
                       type* operator[](int index){
                             return (type*) table[index];
                             }
                             
                       //! Copies row 'r' into pre-allocated row vector 'target'
                       void copyIntoArray(type *target, int r){
                            for(int i=0; i<cols; i++){
                                    target[i]=table[r][i];
                                    }
                            }
                            
                       //! Copies 'table' into pre-allocated 'target'
                       void copyIntoTable(type **target){
                            for(int r=0; r<rows; r++){
                                     for(int c=0; c<cols; c++){
                                              target[r][c]=table[r][c];
                                              }
                                     }
                            }
                            
                       //! Returns pointer to a copy of this with an embedded copy of data
                       matrix *clone(){
                              type **nn=(type **)calloc(rows,sizeof(type *));
                              for(int i=0; i<rows; i++){
                                      nn[i]=(type *)calloc(rows,sizeof(type));
                                      copyIntoArray(nn[i],i);
                                      }
                              matrix *vv=new matrix(nn,rows,cols);
                              vv->embed();
                              }
                              
                       //! Matrix operations
                       void add(matrix *v){
                            for(int r=0; r<min(rows,v->rows); r++){
                                     for(int c=0; c<min(cols,v->cols); c++){
                                              table[r][c]+=(*v)[r][c];
                                              }
                                     }
                            }
                       void sub(matrix *v){
                            for(int r=0; r<min(rows,v->rows); r++){
                                     for(int c=0; c<min(cols,v->cols); c++){
                                              table[r][c]-=(*v)[r][c];
                                              }
                                     }
                            }
                       void inverse(){
                            for(int r=0; r<rows; r++){
                                     for(int c=0; c<cols; c++){
                                              table[r][c]=-table[r][c];
                                              }
                                     }
                            }
                            
                       //! Matrix-vector product: returns pointer to a vector with embedded data
                       vector<type> *prod(vector<type> *v){
                            if(cols!=v->length) return NULL;
                            type *result=(type *)calloc(rows,sizeof(type));
                            int  hh=cols;
                            for(int r=0; r<rows; r++){
                                     result[r]=0;
                                      for(int h=0; h<hh; h++){
                                               result[r]+=table[r][h]*(*v)[h];
                                               }
                                     }
                            vector<type> *vv=new vector<type>(result,rows);
                            vv->embed();
                            return vv;
                            }
                            
                       //! Matrix-vector product: returns pointer to a vector with embedded data
                       vector<type> operator*(vector<type> v){
                            type *result=(type *)calloc(rows,sizeof(type));
                            int  hh=(cols<v.length?cols:v.length);
                            for(int r=0; r<rows; r++){
                                     result[r]=0;
                                      for(int h=0; h<hh; h++){
                                               result[r]+=table[r][h]*(v[h]);
                                               }
                                     }
                            vector<type> vv(result,rows);
                            vv.embed();
                            return vv;
                            }   
                            
                       //! Matrix-matrix product: returns pointer to a matrix with embedded data
                       matrix *prod(matrix *m){
                            if(cols!=m->rows) return NULL;
                            type **result=(type **)calloc(rows,sizeof(type *));
                            for(int r=0; r<rows; r++){
                                    result[r]=(type *)calloc(m->cols,sizeof(type));
                                    }
                            int  hh=cols;
                            int cc=m->cols;
                            
                            for(int r=0; r<rows; r++){
                                    //cerr<<"r: "<<r<<endl;
                                    for(int c=0; c<cc; c++){
                                    //cerr<<"c: "<<c<<endl;
                                     result[r][c]=0;
                                      for(int h=0; h<hh; h++){
                                    //cerr<<"h: "<<h<<endl;
                                               result[r][c]+=table[r][h]*(*m)[h][c];
                                               }
                                      }
                                     }
                            matrix *vv=new matrix(result,rows,cc);
                            vv->embed();
                            return vv;
                            }
                            
                        //! Matrix-matrix product: returns a matrix with embedded data
                       matrix operator*(matrix m){
                            if(cols!=m.rows) cerr<<"Matrici di dimensioni differenti!!!"<<endl;
                            int cols=min(cols,m.rows);
                            type **result=(type **)calloc(rows,sizeof(type *));
                            for(int r=0; r<rows; r++){
                                    result[r]=(type *)calloc(m.cols,sizeof(type));
                                    }
                            int  hh=cols;
                            int cc=m.cols;
                            
                            for(int r=0; r<rows; r++){
                                    //cerr<<"r: "<<r<<endl;
                                    for(int c=0; c<cc; c++){
                                    //cerr<<"c: "<<c<<endl;
                                     result[r][c]=0;
                                      for(int h=0; h<hh; h++){
                                    //cerr<<"h: "<<h<<endl;
                                               result[r][c]+=table[r][h]*(m)[h][c];
                                               }
                                      }
                                     }
                            matrix vv(result,rows,cc);
                            vv.embed();
                            return vv;
                            }
                            
                       //! Prints the matrix contents     
                             void print(){
                                  for(int r=0; r<rows; r++){
                                          for(int c=0; c<cols; c++){
                                                  cerr<<table[r][c]<<"\t";
                                                  }
                                                  cerr<<endl;
                                          }
                                          
                                  }
                                  
                       //! Summation of row 'rr'
                       type rowSum(int rr){
                            type sum=0;
                            for(int c=0; c<cols; c++){
                                    sum+=table[rr][c];
                                    }
                            return sum;
                            }
                       
                       //! Summation of col 'cc'
                       type colSum(int cc){
                            type sum=0;
                            for(int r=0; r<rows; r++){
                                    sum+=table[r][cc];
                                    }
                            return sum;
                            }
                            
                       //! Summation of absolute values of row 'rr' elements
                       type rowAbsSum(int rr){
                            type sum=0;
                            for(int c=0; c<cols; c++){
                                    sum+=(table[rr][c]>=0?table[rr][c]:-table[rr][c]);
                                    }
                            return sum;
                            }
                       
                       
                       //! Summation of absolute values of col 'cc' elements
                       type colAbsSum(int cc){
                            type sum=0;
                            for(int r=0; r<rows; r++){
                                    sum+=(table[r][cc]>=0?table[r][cc]:-table[r][cc]);
                                    }
                            return sum;
                            }
                            
                       //! 1-norm of the matrix
                       type norm_1(){
                            type norm=0;
                            for(int c=0; c<cols; c++){
                                    norm=max(norm,colAbsSum(c));
                                    }
                            return norm;
                            }
                            
                       //! Infinite-norm of the matrix
                       type norm_INF(){
                            type norm=0;
                            for(int r=0; r<rows; r++){
                                    norm=max(norm,rowAbsSum(r));
                                    }
                            return norm;
                            }
                            
                       //! 2-norm of the matrix
                       type norm_frobenius(){
                            type norm=0;
                            for(int r=0; r<rows; r++){
                                    for(int c=0; c<cols; c++){
                                            norm+=table[r][c]*table[r][c];
                                            }
                                    }
                            return sqrt(norm);
                            }
                            
                       //! Summation of all the elements in the matrix
                       type summation(){
                            type res=0;
                            for(int r=0; r<rows; r++){
                                    for(int c=0; c<cols;c++){
                                            res+=table[r][c];
                                            }
                                    }
                            return res;
                            }
                            
                       //! Mean of the matrix
                       type mean(){
                            type res=summation();
                            return res/(rows*cols);
                            }
                };
                
                //! Gauss elimination algorithm
                //! Warning: it works over the given matrix and vector, so the data inside will be changed during the computation
                //! Warning: the result is not guaranteed to be a valid, consistent data array unless the ending status is 'finite'
                template <class type>
                class gauss{
                      private:
                              matrix<type> *m; //! Coefficient matrix
                              vector<type> *b; //! Known terms vector
                              int *indexs; //! Auxiliary vector to support partial pivoting
                              vector<type> *result; //! Result containing vector (undefined in case of non completion)
                              char flag; //! Completion status:
                                   //! f: one, finite solution
                                   //! n: no solution
                                   //! i: infinite solutions
                                   //! k: (during computation) the algorithm is proceding without errors
                                   //! e: there has been errors
                      public:
                             //! Builds a gauss algorithm over the coefficient matrix 'mym' and the known terms vector 'myb'
                             gauss(matrix<type> *mym, vector<type> *myb){
                                                //! Checks the rows of the matrix and the length of the column vector to be the same number
                                                if(mym->rows==myb->length){
                                                m=mym;
                                                b=myb;
                                                flag='k';
                                                result=NULL;
                                                indexs=(int *) calloc(mym->rows,sizeof(int));
                                                for(int c=0; c<m->rows; c++){
                                                        indexs[c]=c;
                                                        }
                                                }
                                                //! If the rows/length check fails, an error condition is signalled
                                                else{
                                                     m=NULL;
                                                     b=NULL;
                                                     indexs=NULL;
                                                     flag='e';
                                                     result=NULL;
                                                     }
                                                }
                             //! Flag setters
                             void setFinite(){
                                  flag='f';
                                  }
                             void setNull(){
                                  flag='n';
                                  }
                             void setInfinite(){
                                  flag='i';
                                  }
                             void setError(){
                                  flag='e';
                                  }
                             void setOK(){
                                  flag='k';
                                  }
                                  
                             //! Status getters
                             bool isFinite(){ //! one, finite solution exists in the result vector
                                  return flag=='f';
                                  }
                             bool isNull(){ //! there is no solution
                                  return flag=='n';
                                  }
                             bool isInfinite(){ //! there are infinite solutions
                                  return flag=='i';
                                  }
                             bool isError(){ //! there has been an error computing
                                  return flag=='e';
                                  }
                             bool isOK(){ //! the algorithm is proceeding well, or it's ended with a single, finite solution
                                  return flag=='k'||isFinite();
                                  }
                             
                             //! Pointer to the solution of the algorithm, or NULL in case of error/non completion
                             vector<type> *solution(){
                                  if(isFinite()) return result;
                                  else return NULL;
                                  }
                                  
                             //! Reference to an element of the complete matrix
                             type& el(int r, int c){
                                  if(r>=0&&r<m->rows&&c>=0&&c<(m->cols+1)){
                                  if(c==m->cols) return (*b)[indexs[r]];
                                  else return (*m)[indexs[r]][c];                                         
                                  }
                                  else setError();
                                  return (*m)[0][0];
                                  }
                                  
                             //! Reference to an element of the known terms vector
                             type& elB(int r){
                                  return el(r,m->cols);
                                  }
                             
                             //! Row swap elementary operation     
                             void swap(int i, int j){
                                  int tmp=indexs[j];
                                  indexs[j]=indexs[i];
                                  indexs[i]=tmp;
                                  }
                             
                             //! Row multiplication elementary operation     
                             void mult(int r, type val){
                                  for(int c=0; c<=m->cols; c++){
                                          el(r,c)*=val;
                                          }
                                  }
                                  
                             //! Constant multiplyed row addition elementary operation
                             void add(int target, int origin, type constant){
                                  if(constant==0) return;
                                  for(int c=0; c<=m->cols; c++){
                                          el(target,c)+=constant*el(origin,c);
                                          }
                                  }
                             
                             //! Value of the pivot in row 'r'
                             type pivotVal(int r){
                                  for(int c=0; c<m->cols; c++){
                                          type elem=el(r,c);
                                          if(elem!=0) return elem;
                                          }
                                  return 0;
                                  }
                             
                             //! Position of the pivot in row 'r'
                             int pivotPos(int r){
                                  for(int c=0; c<m->cols; c++){
                                          type elem=el(r,c);
                                          if(elem!=0) return c;
                                          }
                                  return -1;
                                  }
                                 
                             //! Position of the pivot from right in row 'r'
                             int reversePivotPos(int r){
                                  for(int c=m->cols-1; c>=0; c--){
                                          type elem=el(r,c);
                                          if(elem!=0) return c;
                                          }
                                  return -1;
                                  }
                                  
                             //! Uses 'add' operation to set to 0 all the elements in column 'col' under the one in row 'beginningRow', 
                             //!      which has pivot value 'pivot'; 
                             //!      then returns the number of the row with maximum absolute vlued pivot for partial pivoting.
                             int eliminateCol(int beginningRow, int col, type pivot){
                                 type max=0;
                                 int maxPivotRow=beginningRow;
                                 //! Eliminates all the elements in the given column...
                                 for(int r=beginningRow+1; r<m->rows; r++){
                                         add(r,beginningRow,-el(r,col)/pivot);
                                         el(r,col)=0;
                                         if(col+1<m->cols){
                                         type tmp=el(r,col+1);
                                         type tmpa=(tmp<0 ? -tmp : tmp);
                                         if(tmpa>max){
                                                      max=tmpa;
                                                      maxPivotRow=r;
                                                      }}
                                         }
                                 int cc=col+1;
                                 //! ...then looks for the row with the maximum absolute valued pivot element and returns it
                                 while(maxPivotRow==beginningRow){
                                                                  for(int r=beginningRow+1;r<m->rows;r++){
                                                                          if(cc+1<m->cols){
                                                                          type tmp=el(r,cc+1);
                                                                          type tmpa=(tmp<0 ? -tmp : tmp);
                                                                          if(tmpa>max){
                                                                                       max=tmpa;
                                                                                       maxPivotRow=r;
                                                                                       }}
                                                                          }
                                                                          if(cc>=m->cols) break;
                                                                          cc++;
                                                                  }
                                 return maxPivotRow;
                                 }
                                 
                             //! Transforms the complete matrix in echelon form, and returns true if the operation succedes
                             bool toEchelon(){
                                  if(!isOK())return false;
                                  for(int r=0; r<m->rows; r++){ //! For every row
                                          int pivCol=pivotPos(r); //! Finds pivot position
                                          if(pivCol==-1&&elB(r)!=0){ //! Aborts in case of 0=1
                                                                    setNull();
                                                                    return false;
                                                                    }
                                          if(pivCol==-1) continue; //! If no pivot could be found, cicles until the end of matrix: further passages will set eventually set the 'infinite' status
                                          
                                          type pivVal=el(r,pivCol); //! Saves the pivot value
                                          int maxPRow=eliminateCol(r,pivCol,pivVal); //! Nullifies the values under the pivot
                                          if(r+2<m->rows&&maxPRow!=r) swap(r+1,maxPRow); //! Partial pivoting
                                          }
                                  return true;
                                  }
                                  
                             //! Given the complete matrix in echelon form, it diagonalizes the coefficient matrix and returns true if the operation succedes.
                             bool diagonalize(){
                                  int cc=m->cols-1;
                                  for(int r=m->rows-1; r>=0; r--){
                                          if(!isOK()){
                                                        return false;
                                                        }
                                          if(r<cc){ //! If there is a null columns, the status is set to 'infinite'
                                                   setInfinite();
                                                   return false;
                                                   }
                                          int mypivotcol=reversePivotPos(r); //! Finds the pivot from right
                                          type mypivot=el(r,mypivotcol); //! Saves the pivot value
                                          for(int rr=r-1; rr>=0; rr--){ //! Nullifies all the values upwards
                                                  add(rr,r,-el(rr,mypivotcol)/mypivot);
                                                  el(rr,mypivotcol)=0;
                                                  }
                                          cc--;
                                          }
                                  return true;
                                  }
                                  
                             //! Given the diagonalized coefficient matix, it fills the result vector, and returns true if the operation succedes
                             bool fillVector(){
                                  if(!isOK()) return false;
                                  for(int r=0; r<m->rows; r++){
                                          int mypivotcol=pivotPos(r);
                                          if(mypivotcol<0){ //! if there is a null row, sets the status to 'infinite'
                                                           setInfinite();
                                                           return false;
                                                           }
                                          type mypivot=el(r,mypivotcol);
                                          (*result)[r]=elB(r)/mypivot;
                                          }
                                  return true;
                                  }
                             
                             //! Product of the diagonal elements
                             type diagonalProd(){
                                  type res=1;
                                  for(int r=0; r<m->rows; r++){
                                            res*=el(r,r);
                                            }
                                  return res;
                                  }
                             
                             //! Returns the determinant of the diagonalized form
                             type determinant(){
                                  if(isFinite()){
                                                 return diagonalProd();
                                                 }
                                  else if(isInfinite()||isNull()){
                                       return 0;
                                       }
                                  }
                             
                             //! Prints the matrix contents     
                             void print(){
                                  if(isOK()){
                                  for(int r=0; r<m->rows; r++){
                                          for(int c=0; c<=m->cols; c++){
                                                  cout<<el(r,c)<<"\t";
                                                  }
                                                  cout<<endl;
                                          }}
                                          cout<<"Err:"<<isError()<<"Fin:"<<isFinite()<<"Nul:"<<isNull()<<"Inf:"<<isInfinite()<<endl;
                                          if(isFinite()){cout<<"Result: "; result->print();}
                                          //cout<<endl;
                                  }
                             
                             //! Executes the gaussian elimination algorithm
                             void solve(){
                                  //! Matrix to echelon form
                                          if(!isOK())return;
                                          int r=toEchelon();
                                  //! Matrix diagonalization
                                          if(!isOK())return;
                                          r+=diagonalize();
                                  //! Result vector filling and 'finite' status setting
                                          if(!isOK())return;
                                          result=new vector<type>(b->length);
                                          r+=fillVector();
                                          if(r==3) setFinite(); //! if r==3, it means that the three operations were all succesfull, and the algorithm is thus completed
                                          
                                  }
                      };
                
          }
