# include "bucket.hpp"

int Bucket::N = 0;

inline int inc(int dir, int incr){
    if (incr < 0) throw("invalid increment in inc, must be postive");
    return (dir+incr)%8;
}

inline int Bucket::dist(const shared_ptr<Bucket> other) const {
    return (ind_i-other->i())*(ind_i-other->i()) + 
        (ind_j-other->j())*(ind_j-other->j());
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
void Bucket::addToCorner(int diag){ // to be called from corner!
    if (!(diag%2)) throw ("diag must be uneven in addToCorner");
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


/*
to_go = {cnt dir 0, cnt dir 2, cnt dir 4, cnt dir 6, last step, original dir}
all dirs given from 0 to 7
*/
shared_ptr<Bucket> Bucket::searchCorner(vector<int>& to_go){

    cout << "in search corner " << endl;
    shared_ptr<Bucket> current = self;
    bool cond = current->isCorner() == -1;
    int next_step = -1;
    int possible_step = -1;

    while(cond) { // while current is not a corner
        cout << "in while" << endl;

        next_step = -1;
        for (int k=0; k<4; k++) {
            cout << "in for" << endl;
            /* don't allow to go back*/
            if (k==inc(to_go[4],4)/2) continue; 
            /* don't allow to go "behind" original bucket */
            if (k==inc(to_go[5],4)/2 && to_go[k]<0) continue;

            if (to_go[k] > 0 && current->neighbours[k*2] != nullptr) {
                next_step = k*2;
                break;
            }
            possible_step = 2*k;

        }

        // cout << "after for " << endl;

        if (next_step == -1)
            next_step = possible_step;

        to_go[4] = next_step; // update last step
        to_go[next_step/2]--;
        to_go[inc(next_step, 4)]++;
        current = current->neighbours[next_step];

        cout << "cond: " << (current->isCorner() == -1) << endl;
        cond = current->isCorner() == -1;
    }

    return current;
}

/* test function to be deleted in the end */
void Bucket::test() { // assume N = 10, (i,j) = (5,5)
    shared_ptr<Bucket> current = self;
    current->newBucket(6);
    while (current->neighbours[6] != nullptr) {
        current = current->neighbours[6];
        current->newBucket(6);
    }

    current = neighbours[6];
    current = current->neighbours[6];
    current = current->neighbours[0];

    current->newBucket(0);
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
        if (dir%2) addToCorner(dir);
        else addToCorner(inc(dir,1));
    } else if (diag == -1) { // if not a corner

        cout << "not a corner" << endl;
        // throw("not a corner in newBucket, still to be implemented"); 
        switch (dir)
        {
        case 1:
            neighbours[0]->newBucket(2);
            return;
            break;
        case 3:
            neighbours[4]->newBucket(2);
            return;
            break;
        case 5:
            neighbours[4]->newBucket(6);
            return;
            break;
        case 7:
            neighbours[0]->newBucket(6);
            return;
            break;
        default:
            break;
        }

        if (dir%2) throw("dir must be even if not a corner");

        /* search nearest corner and write it in 'corner' 
           note that accessed neighbours cannot be nullptr
        */
        vector<int> to_go(6,0);
        shared_ptr<Bucket> corner;
        to_go[dir/2]++;
        to_go[inc(dir,2)/2]++;
        to_go[inc(dir,4)/2]--;
        to_go[inc(dir,6)/2]--; 
        to_go[4] = inc(dir,2); // last step
        to_go[5] = dir; // original dir to add 
        corner = neighbours[inc(dir,6)]->searchCorner(to_go); 
        {
            vector<int> to_go_tmp(6,0);
            shared_ptr<Bucket> tmp;
            to_go_tmp[dir/2]++;
            to_go_tmp[inc(dir,6)/2]++; 
            to_go_tmp[inc(dir,4)/2]--;
            to_go_tmp[inc(dir,2)/2]--; 
            to_go_tmp[4] = inc(dir,6); // last step
            to_go_tmp[5] = dir; // original dir to add 
            tmp = neighbours[inc(dir,2)]->searchCorner(to_go_tmp); 
            if (dist(tmp) < dist(corner)) {
                corner = tmp;
                to_go = move(to_go_tmp);
            }
        }

        cout << "found corner ("<<corner->ind_i<<","<<corner->ind_j<<")" << endl;
 
    } else {
        addToCorner(diag);
    }


    // /* just for testing */
    // if (neighbours[dir] != nullptr) {
    //     neighbours[dir]->newBucket(6); 
    // } // WORKING :DDD
}