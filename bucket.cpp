# include "bucket.hpp"

int Bucket::N = 0;

void Bucket::newBucket(int dir){
    if (dir<0 || dir>7) throw ("invalid direction in newBucket");
    /* compute global indices for new bucket */
    int i,j;
    if (0<dir && dir<4) i = ind_i+1;
    if (4<dir) i = ind_i-1;
    if (2<dir && dir<6) j = ind_j-1;
    if (2>dir || dir>6) j = ind_j+1;

    if (i==N || j==N || i<0 || j<0) {
        cout << "Index out of bound in newBucket. No new Bucket created" << endl;
        return;
    }
    if (neighbours[dir] != nullptr) {
        cout << "Neighbour already exists" << endl;
        return;
    }

    neighbours[dir] = new Bucket(i,j);
    neighbours[dir]->neighbours[inc(dir, 4)] = this;

    bool go_cw = 0;
    int cd = dir; // current direction of new Bucket (viewpoint current bucket)
    int next_dir; // dir to go 
    // int assign; // assign direction viewpoint current bucket
    Bucket* current = this;
    for (int k = 1; k<8; k++){

        /* go in previous even dir */
        if(inc(cd,7)%2) next_dir = inc(cd, 6);
        else next_dir = inc(cd, 7);

        /* check if neighbour next_dir exists */
        if (current->neighbours[next_dir] == nullptr) {
            go_cw = 1;

        }
        /* update current dir of bucket */
        cd = inc(cd,1);
        

        


    }
    
}