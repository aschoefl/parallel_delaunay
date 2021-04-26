# include "bucket.hpp"

shared_ptr<Bucket> Bucket::bb = move(shared_ptr<Bucket>((Bucket*) new BoundaryBucket()));
list<shared_ptr<Bucket>> Bucket::buckets;
shared_ptr<Bucket> Bucket::root = nullptr;
int Bucket::P=0 , Bucket::N=0;
double Bucket::buffer[MAX_PNTS*2+1];


void Bucket::addToList(){
    if (!(find(buckets.begin(), buckets.end(), self) != buckets.end()))
        buckets.push_back(self);
    // Point p(ind_i, ind_j);
    // if (!(find(buckets.begin(), buckets.end(), p) != buckets.end()))
    //     buckets.push_back(move(p));
}

void Bucket::printList(){
    ofstream myfile ("points"+to_string(r())+".txt");
    if (myfile.is_open()) {
        myfile << "[";
        for (auto b : buckets) {
            for (auto p: b->points)
                myfile << p << ", ";
        }
        myfile << "]" << endl;
        myfile.close();
    }
    else cout << "Unable to open file";

    // buckets.sort();
    // cout << "List of buckets [";
    // for (auto p : buckets) 
    //     cout << p << ", ";
    // cout << "]" << endl;

}

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

    if (dir<0 || dir>7) throw ("invalid direction in addBucket");
    if (neighbours[dir] != nullptr) {
        cout << "Neighbour already exists" << endl;
        return;
    }

    /* compute global indices for new bucket */
    int i,j;
    getIndex(dir, i,j);

    if (i==N || j==N || i<0 || j<0) {
        cout << "Indices ("<<i<<","<<j<<") out of bound in addBuket. No Bucket added" << endl;
        return;
    }
    // cout << "add Bucket with global indices ("<<i<<","<<j<<") in dir " << dir << endl;
    
    int cd = dir; // current direction of new Bucket (viewpoint current bucket)
    int next_dir; // dir to go 
    shared_ptr<Bucket> current = self;
    /* work around the memory leaks */
    shared_ptr<Bucket> new_bucket = (new Bucket(i,j))->self;

    // just for testing, but not working
    // new_bucket->fillCoordinates();

    for (auto k = 0; k<8; k++){ 
        /* assign boundary buckets */
        new_bucket->getIndex(k,i,j);
        if (i==N || j==N || i<0 || j<0)
            new_bucket->neighbours[k] = bb;
    }

    for (auto k=0; k<8; k++) {

        /* set links to each other */
        current->neighbours[cd] = new_bucket; 
        new_bucket->neighbours[inc(cd, 4)] = current->self;

        /* go in previous even dir */
        if(inc(cd,7)%2) next_dir = inc(cd, 6);
        else next_dir = inc(cd, 7);
        
        if (current->neighbours[next_dir] == nullptr || current->neighbours[next_dir]->isBnd()) {
            if (inc(cd,7)%2) { // if previous dir is uneven
                next_dir = inc(cd, 7);
                if (current->neighbours[next_dir] == nullptr || current->neighbours[next_dir]->isBnd())
                    break;
                /* skip one iteration */
                k++;
                cd = inc(cd,1); 
            } else break;
        } 
        
        current = current->neighbours[next_dir];
        cd = inc(cd,1);
    }
    cd = dir; 
    current = self;
    for (int k = 1; k<8; k++){ //go clockwise until nullptr
        /* go in next even dir */
        if(inc(cd,1)%2) next_dir = inc(cd, 2);
        else next_dir = inc(cd, 1);

        /* check if neighbour next_dir exists */
        if (current->neighbours[next_dir] == nullptr || current->neighbours[next_dir]->isBnd()) {
            if (inc(cd,1)%2) {
                next_dir = inc(cd,1);
                if (current->neighbours[next_dir] == nullptr || current->neighbours[next_dir]->isBnd())
                    break;
                /* skip one iteration */
                k++;
                cd = inc(cd,7); 
            } else break;
        }
        current = current->neighbours[next_dir];

        /* update dir of new bucket */
        cd = inc(cd,7);

        /* set links to each other */
        current->neighbours[cd] = new_bucket; 
        new_bucket->neighbours[inc(cd, 4)] = current->self;
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
                if (self == root) {
                    current->addBucket(k);
                    next = current->neighbours[k];
                } else {
                    return (*root)(i,j);
                }
            }
            if (next->isBnd())
                throw runtime_error("index out of bound in ()");
            current = next;
            dir[k]--;
        }
    }
    return current;
}

void Bucket::addPoint(double x, double y) {
    Point pnt(x,y);
    auto pos = find_if(points.begin(), points.end(), [pnt](auto s) {
        return s < pnt;
    });
    points.insert(pos, move(pnt));
}

void Bucket::fillCoordinates(){

    if (!(coordinates.empty())) return;

    if (r() != root->r()) {
        /* get data from right processor */
        int no_bucket = i()*N+j();
        MPI_Send(&no_bucket,1,MPI_INT,r(),TAG_NOTIFY,MPI_COMM_WORLD);
        MPI_Recv(buffer,MAX_PNTS,MPI_DOUBLE,r(),TAG_DATA,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        // MPI_Request req;
        // int no_bucket = i()*N+j();
        // MPI_Isend(&no_bucket,1,MPI_INT,r(),TAG_NOTIFY,MPI_COMM_WORLD, &req);
        // MPI_Recv(buffer,MAX_PNTS,MPI_DOUBLE,r(),TAG_DATA,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        // MPI_Request_free(&req);

        coordinates.push_back(buffer[0]);
        for (auto k=0; k<2*buffer[0]; k+=2) {
            coordinates.push_back(buffer[k+1]);
            coordinates.push_back(buffer[k+2]);
            addPoint(buffer[k+1],buffer[k+2]);
        }
        // throw runtime_error("get from other processer not yet implemented");
        
        cout << "filled coordinates of (" << i()<<", "<<j()<<
        ") on processor "<<r()<<endl;
        return;
    }

    /* number of points in this bucket 
       pfusch distribution  */
    int cnt = rand() % (MAX_PNTS+1);
    for (auto i=0; i<10; i++)
        cnt = max(min(cnt, rand() % (MAX_PNTS+1)),1);
    
    auto xmin = i()/static_cast<double>(N);
    auto xmax = (i()+1)/static_cast<double>(N);
    auto ymin = j()/static_cast<double>(N);
    auto ymax = (j()+1)/static_cast<double>(N);
    coordinates.push_back(cnt); // set amount of points 

    /* assign random coordinates in right bucket */
    for (auto i=0; i<cnt; i++){
        auto x = xmin + static_cast<double>(rand())/( static_cast<double>(RAND_MAX/(xmax-xmin)));
        auto y = ymin + static_cast<double>(rand())/( static_cast<double>(RAND_MAX/(ymax-ymin)));
        addPoint(x,y);
        coordinates.push_back(x);
        coordinates.push_back(y);
    }

    // cout << "filled coordinates of (" << i()<<", "<<j()<<
    // ") on processor "<<r()<<endl;

}


void Bucket::sendCoordinates(int destination, int no_bucket) {

    auto ii = no_bucket/N;
    auto jj = no_bucket%N; 

    if (i() != ii || j() != jj ) {
        (*root)(ii,jj)->sendCoordinates(destination, no_bucket);
        return;
    }

    /* get or gerenate data if not available */
    if (coordinates.empty())
        fillCoordinates();

    MPI_Send(coordinates.data(), MAX_PNTS, MPI_DOUBLE, destination, 
                TAG_DATA, MPI_COMM_WORLD);
}

vector<Point> Bucket::getPoints() {
    /* get or gerenate data if not available, 
    check if coordinates are empty because points can
    remain empty if no point in bucket */
    if (coordinates.empty())
        fillCoordinates();
    return points;
}



/* source for how to provide fairness:
https://www.mcs.anl.gov/research/projects/mpi/tutorial/gropp/node93.html#Node93
signatures

MPI_Send(
    void* data,
    int count,
    MPI_Datatype datatype,
    int destination,
    int tag,
    MPI_Comm communicator)

MPI_Recv(
    void* data,
    int count,
    MPI_Datatype datatype,
    int source,
    int tag,
    MPI_Comm communicator,
    MPI_Status* status)

*/
void Bucket::doSomething() {
    if (running) throw runtime_error("already running");
    running = true;

    MPI_Request requests[MAX_PROC]; 
    MPI_Status  statuses[MAX_PROC]; 
    int         indices[MAX_PROC]; 
    int         buf[MAX_PROC]; 
    auto size = P*P;
    int ndone;
    for (auto i=0; i<size; i++)  
        MPI_Irecv( buf+i, 1, MPI_INT, i, 
                TAG_NOTIFY, MPI_COMM_WORLD, &requests[i] ); 
    bool startup = true;
    while(running) { 
        if (startup) {
            test();
            startup = false;
        }
        MPI_Waitsome( size, requests, &ndone, indices, statuses ); 
        for (auto i=0; i<ndone; i++) { 
            auto j = indices[i]; 
            printf( "Msg from %d with tag %d and value %d\n",  
                    statuses[i].MPI_SOURCE,  
                    statuses[i].MPI_TAG,
                    buf[j]); 
            if (statuses[i].MPI_SOURCE != r())
                sendCoordinates(statuses[i].MPI_SOURCE, buf[j]);
            // else { /* stop if receive msg from oneself */
            //     cout << "stop running " << endl;
            //     running = false;
            // }
            MPI_Irecv( buf+j, 1, MPI_INT, j, 
                    TAG_NOTIFY, MPI_COMM_WORLD, &requests[j]); 
        } 
    } 
}
void Bucket::calculateDelauney(){

}
/* test function to be deleted in the end */
void Bucket::test() { 

    shared_ptr<Bucket> a = self;
    a = (*a)(5,5);
    // cout << "here " << endl;
    // a = (*a)(min(i()+5, N-1),min(i()+5, N-1));  // this is a problem
    // a = (*a)(N-1,N-1); 
    // a = (*a)(5,5); 
    // a = (*a)(5,5);
    // a = (*a)(0,7);
    // a = (*a)(1,N-2);
    // a = (*a)(3,1);
    // a = (*a)(7,7);

    // problem: deadlock somewhere
    for (int ii =5;ii<N;ii++) {
        for(int jj=5;jj<N;jj++) {
            cout << r() << "in for" << endl;
            a = (*a)(ii,jj);
            a->getPoints();
        }
    }

    // cout << "r: " << a->r() << endl;
    vector<Point> pnts = move(a->getPoints());
    // cout << "points on " << r() << " ";
    // for (auto p: pnts)
    //     cout << p << " ";
    // cout << endl;

    // int no_bucket = i()*N+j();
    // MPI_Send(&no_bucket,1,MPI_INT,r(),TAG_NOTIFY,MPI_COMM_WORLD);


    printList();
    cout << "****** Proc " << r() << " finished ******" << endl;
    // cout << "coordinates: ";
    // for (auto p: coordinates)
    //     cout << p << " ";
    // cout << endl;

    // a->printNeighbours();
    // a->printList();

    /* send stop notification */
    // int stop = -1;
    // MPI_Send(&stop,1,MPI_INT,r(),TAG_NOTIFY,MPI_COMM_WORLD);
    // cout << "here " << r() << endl;
}

