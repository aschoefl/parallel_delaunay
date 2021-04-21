# include "bucket.hpp"

int Bucket::N = 0;

inline int inc(int dir, int incr){
    if (incr < 0) throw("invalid increment in inc, must be postive");
    return (dir+incr)%8;
}

/*  
if bucket is an ordernary corner return free diagonal dir, 
if bucket is a sinle element (not connected to anything) return 8,
else not return -1
*/
int Bucket::isCorner() {
    int cnt = 0;
    int dir[2];
    for (int i=0; i<8; i+=2){ // iterate through even neighbours
        if (neighbours[i] == nullptr) {
            if (cnt == 0) dir[0]=i;
            if (cnt == 1) dir[1]=i;
            cnt++;
        }
    }
    if (cnt==2) { // assign diagonal dir
        if (inc(dir[0],2)==dir[1]) return inc(dir[0],1);
        else return inc(dir[1],1);
    }
    if (cnt==4) return 8;
    return -1;
}

/* diag .. dir in diag direction (viewpoint corner) */
void Bucket::newToCorner(int diag){ // to be called from corner!
    if (!(diag%2)) throw ("diag must be uneven in newToCorner");
    addBucket(inc(diag,1)); 
    addBucket(inc(diag,7));
    addBucket(diag);

    cout << "neighbours of corner ("<<ind_i<<","<<ind_j<<"): ";
    for (auto bucket : neighbours){
        if (bucket != nullptr)
            cout << " ("<<bucket->ind_i<<","<<bucket->ind_j<<") ";
    }
    cout << endl;

    // cout << inc(diag,1) << inc(diag,7) << diag << endl;
}

void Bucket::addBucket(int dir){ // to be called from corner

    /* compute global indices for new bucket */
    int i = ind_i;
    int j = ind_j;
    if (0<dir && dir<4) j++;
    if (4<dir) j--;
    if (2<dir && dir<6) i--;
    if (2>dir || dir>6) i++;

    if (i==N || j==N || i<0 || j<0) {
        cout << "Indices ("<<i<<","<<j<<") out of bound in newBucket. No new Bucket created" << endl;
        return;
    }
    // cout << "add Bucket with global indices ("<<i<<","<<j<<") in dir " << dir << endl;
    
    int cd = dir; // current direction of new Bucket (viewpoint current bucket)
    int next_dir; // dir to go 
    shared_ptr<Bucket> current = self;
    /* work around the memory leaks */
    shared_ptr<Bucket> new_bucket = (new Bucket(i,j))->self;

    for (int k = 0; k<8; k++){ //go counter clockwise until nullptr
        /* set links to each other */
        current->neighbours[cd] = new_bucket; 
        new_bucket->neighbours[inc(cd, 4)] = current->self;

        /* go in previous even dir */
        if(inc(cd,7)%2) next_dir = inc(cd, 6);
        else next_dir = inc(cd, 7);

        /* check if neighbour next_dir exists */
        if (current->neighbours[next_dir] == nullptr) break;
        else current = current->neighbours[next_dir];

        /* update dir of new bucket */
        cd = inc(cd,1);
    }
    cd = dir; 
    current = self;
    for (int k = 1; k<8; k++){ //go clockwise until nullptr
        /* go in next even dir */
        if(inc(cd,7)%2) next_dir = inc(cd, 2);
        else next_dir = inc(cd, 1);

        /* check if neighbour next_dir exists */
        if (current->neighbours[next_dir] == nullptr) break;
        else current = current->neighbours[next_dir];

        /* update dir of new bucket */
        cd = inc(cd,7);

        /* set links to each other */
        current->neighbours[cd] = new_bucket; 
        new_bucket->neighbours[inc(cd, 4)] = current->self;
    }
}


void Bucket::newBucket(int dir){ // public function

    if (dir<0 || dir>7) throw ("invalid direction in newBucket");
    if (neighbours[dir] != nullptr) {
        cout << "Neighbour already exists" << endl;
        return;
    }

    auto diag = isCorner();
    if (diag == 8) { // if root element, check if dir is diagonal
        if (dir%2) newToCorner(dir);
        else newToCorner(inc(dir,1));
    } else if (diag == -1) { // if not a corner
        // TODO: DEAL WITH THAT CASE!
        throw("not a corner in newBucket, still to be implemented");      
    } else {
        newToCorner(diag);
    }

    // /* just for testing */
    // if (neighbours[dir] != nullptr) {
    //     neighbours[dir]->newBucket(6); 
    // } // WORKING :DDD
}