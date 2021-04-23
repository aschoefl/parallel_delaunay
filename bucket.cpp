# include "bucket.hpp"

int Bucket::N = 0;
shared_ptr<Bucket> Bucket::bb = move(shared_ptr<Bucket>((Bucket*) new BoundaryBucket()));

inline int inc(int dir, int incr){
    if (incr < 0) throw runtime_error("invalid increment in inc, must be postive");
    return (dir+incr)%8;
}

inline int Bucket::dist(const shared_ptr<Bucket> other) const {
    return (ind_i-other->i())*(ind_i-other->i()) + 
        (ind_j-other->j())*(ind_j-other->j());
}

inline void Bucket::printNeighbours() const {
    cout << "neighbours of bucket ("<<ind_i<<","<<ind_j<<"): ";
    for (auto bucket : neighbours){
        if (bucket != nullptr) {
            if (bucket->isBnd()) continue;
            cout << " ("<<bucket->ind_i<<","<<bucket->ind_j<<") ";
        }
    }
    cout << endl;
}

/*  
if bucket is an ordernary corner return free diagonal dir, 
if bucket is a sinle element (not connected to anything) return 8,
else not return -1
*/
int Bucket::isCorner() {

    int np_cnt = 0, bnd_cnt = 0, bnd_dir;
    int dir[2];
    for (int i=0; i<8; i+=2){ // iterate through even neighbours
        if (neighbours[i] == nullptr) {
            if (np_cnt == 0) dir[0]=i;
            if (np_cnt == 1) dir[1]=i;
            np_cnt++;
        }
        else if (neighbours[i]->isBnd()){
            if (bnd_cnt == 0) bnd_dir=i;
            bnd_cnt++;
        }
    }

    if (np_cnt==2) { // assign diagonal dir
        if (inc(dir[0],2)==dir[1]) return inc(dir[0],1);
        else return inc(dir[1],1);
    }
    if (np_cnt==1 && bnd_cnt == 1){
        if (inc(dir[0],2)==bnd_dir) return inc(dir[0],1);
        else return inc(bnd_dir,1);
    }

    if (np_cnt==4) return 8;
    return -1;
}

/* diag .. dir in diag direction (viewpoint corner) */
void Bucket::addToCorner(int diag){ // to be called from corner!
    if (!(diag%2)) throw ("diag must be uneven in addToCorner");
    addBucket(inc(diag,1)); 
    addBucket(inc(diag,7));
    addBucket(diag);

    cout << "addToCorner: ";
    printNeighbours();
}


/*
to_go = {cnt dir 0, cnt dir 2, cnt dir 4, cnt dir 6, last step, original dir}
all dirs given from 0 to 7
*/
shared_ptr<Bucket> Bucket::searchCorner(vector<int>& to_go){

    // cout << "in search corner with to_go:";
    // for (auto v : to_go)
    //    cout <<  v << " ";
    // cout << endl;

    shared_ptr<Bucket> current = self;
    bool cond = current->isCorner() == -1;


    int cnt = 0; // just for testing
    while(cond) { // while current is not a corner

        int next_step = -1;
        int possible_step = -1;
        for (int k=0; k<4; k++) {

            /* don't allow to go back*/
            if (k==inc(to_go[4],4)/2) continue; 

            /* don't allow to go "behind" original bucket */
            if (k==inc(to_go[5],4)/2 && to_go[k]<1) continue;

            /* desirable direction to go */
            if (to_go[k] > 0 && current->neighbours[k*2] != nullptr && !(current->neighbours[k*2]->isBnd()) ) {
                next_step = k*2;
                break;
            }

            /* possible direction to go */
            if (current->neighbours[k*2] != nullptr && !(current->neighbours[k*2]->isBnd())) 
                possible_step = 2*k;

        }

        if (possible_step == -1 && next_step == -1) 
            throw runtime_error("nowhere to go in searchCorner");

        /* go in possible direction if desirable dir not available */
        if (next_step == -1) 
            next_step = possible_step;
        
        to_go[4] = next_step; // update last step
        to_go[next_step/2]--;  // update step counter
        to_go[inc(next_step, 4)/2]++; // update step counter
        current = current->neighbours[next_step]; // update current node

        // current->printNeighbours();
        cond = current->isCorner() == -1; // check if corner reached
        cnt++;

        if (cnt == 100) {
            throw runtime_error("no end after 100 search steps");
        }
    }

    cout << "found corner ("<<current->ind_i<<","<<current->ind_j<<")";
    cout <<  " after " << cnt << " steps in searchCorner" << endl;
    return current;
}

void Bucket::getIndex(int const dir, int& i, int& j){
    /* compute global indices for new bucket */
    i = ind_i;
    j = ind_j;
    if (0<dir && dir<4) j++;
    if (4<dir) j--;
    if (2<dir && dir<6) i--;
    if (2>dir || dir>6) i++;
}


void Bucket::addBucket(int dir){ // to be called from corner

    /* compute global indices for new bucket */
    int i,j;
    getIndex(dir, i,j);

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
    for (int k = 0; k<8; k++){ 
        /* assign boundary buckets */
        new_bucket->getIndex(k,i,j);
        if (i==N || j==N || i<0 || j<0)
            new_bucket->neighbours[k] = bb;
    }

    for (int k = 0; k<8; k++){ //go counter clockwise until nullptr

        /* set links to each other */
        current->neighbours[cd] = new_bucket; 
        new_bucket->neighbours[inc(cd, 4)] = current->self;

        /* go in previous even dir */
        if(inc(cd,7)%2) next_dir = inc(cd, 6);
        else next_dir = inc(cd, 7);

        /* check if neighbour next_dir exists */
        if (current->neighbours[next_dir] == nullptr || current->neighbours[next_dir]->isBnd())
            break;
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
        if (current->neighbours[next_dir] == nullptr || current->neighbours[next_dir]->isBnd())
            break;
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
        cout << "is root " << endl;
        if (dir%2) addToCorner(dir);
        else addToCorner(inc(dir,1));
    } else if (diag == -1) { // if not a corner

        cout << "not a corner" << endl;
        // throw("not a corner in newBucket, still to be implemented"); 
        switch (dir)
        {
        case 1:
            if (neighbours[0] == nullptr) {
                neighbours[2]->newBucket(0);
                return;
            }
            neighbours[0]->newBucket(2);
            return;
            break;
        case 3:
            if (neighbours[4] == nullptr)  {
                neighbours[2]->newBucket(4);
                return;
            }
            neighbours[4]->newBucket(2);
            return;
            break;
        case 5:
            if (neighbours[4] == nullptr)  {
                neighbours[6]->newBucket(4);
                return;
            }
            neighbours[4]->newBucket(6);
            return;
            break;
        case 7:
            if (neighbours[0] == nullptr)  {
                neighbours[6]->newBucket(0);
                return;
            }
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

        int i,j; 
        getIndex(dir, i, j);
        
        auto d = (i-corner->i())*(i-corner->i()) + (j-corner->j())*(j-corner->j());
        while (d > 1) {
            corner->addToCorner(corner->isCorner());
            corner = corner->searchCorner(to_go);
            d = (i-corner->i())*(i-corner->i()) + (j-corner->j())*(j-corner->j());
        }
        corner->addToCorner(corner->isCorner());
 
    } else {
        addToCorner(diag);
    }
}

shared_ptr<Bucket> Bucket::operator() (int i, int j) const{
    auto x = i-ind_i;
    auto y = j-ind_j;
    vector<int> dir(8,0);

    // cout << "x: " << x << " y: " << y << endl;

    /* write as directions */
    if (x>=0 && y>=0) {
        dir[1] = min(x,y);
        dir[0] = x-dir[1];
        dir[2] = y-dir[1];
    }
    if (x<=0 && y>=0) {
        dir[3] = min(-x,y);
        dir[4] = -(x+dir[3]);
        dir[2] = y-dir[3];
    }
    if (x<=0 && y<=0) {
        dir[5] = min(-x,-y);
        dir[4] = -(x+dir[5]);
        dir[6] = -(y+dir[5]);
    }
    if (x>=0 && y<=0) {
        dir[7] = min(x,-y);
        dir[0] = x-dir[7];
        dir[6] = -(y+dir[7]);
    }

    // cout << "dir: ";
    // for (auto v : dir)
    //    cout <<  v << " ";
    // cout << endl;

    /* go to bucket of index and add if needed */
    shared_ptr<Bucket> current = self;
    shared_ptr<Bucket> next;
    for(int k=0; k<8; k++) {
        while (dir[k]) {
            next = current->neighbours[k];
            if (next == nullptr) {
                current->newBucket(k);
                next = current->neighbours[k];
            }
            else if (next->isBnd())
                throw runtime_error("index out of bound in ()");
            current = next;
            dir[k]--;
        }
    }
    return current;
}


/* test function to be deleted in the end */
void Bucket::test() { // assume N = 10, (i,j) = (5,5)
    shared_ptr<Bucket> current = self;
    current->newBucket(6);
    while (current->neighbours[6] != nullptr) {
        if (current->neighbours[6]->isBnd()) break;
        current = current->neighbours[6];
        current->newBucket(6);
        current->neighbours[6]->isBnd();
    }

    shared_ptr<Bucket> a = (*self)(9,8);
    a->printNeighbours();

    auto b = (*a)(0,0);
    b->printNeighbours();

    // current = neighbours[6];
    // current = current->neighbours[6];
    // current = current->neighbours[6];
    // current = current->neighbours[0];

    // current->newBucket(1);
}

