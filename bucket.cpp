# include "bucket.hpp"

shared_ptr<Bucket> Bucket::bb = move(shared_ptr<Bucket>((Bucket*) new BoundaryBucket()));
list<shared_ptr<Bucket>> Bucket::buckets;
shared_ptr<Bucket> Bucket::root = nullptr;
int Bucket::P=0 , Bucket::N=0;
double Bucket::buffer[MAX_PNTS*2+1];

bool evenstep = false;


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

// void Bucket::getIndex(Point const p, int& i, int& j){
//     i = p.x*N;
//     j = p.y*N;
// }

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

int Bucket::fillCoordinates(int stat){

    // MPI_Request req;

    if (!(coordinates.empty())) return 0;

    if (true && r() != root->r()) {

        if (stat == 0) {
            /* get data from right processor */
            int no_bucket = i()*N+j();
            // cout << "send from " << root->r() << " to " << r() << endl;
            MPI_Send(&no_bucket,1,MPI_INT,r(),TAG_NOTIFY,MPI_COMM_WORLD);//, &request);
            return 1;
        } else if (stat == 1) {
            // goto request;
            // resume:
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
            
            // cout << "filled coordinates of (" << i()<<", "<<j()<<
            // ") on processor "<<r()<<endl;
            return 0;
        } else throw runtime_error("invalid status in fillCoordinates");
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
    // for (auto i=0; i<cnt; i++){
    //     auto x = xmin + static_cast<double>(rand())/( static_cast<double>(RAND_MAX/(xmax-xmin)));
    //     auto y = ymin + static_cast<double>(rand())/( static_cast<double>(RAND_MAX/(ymax-ymin)));
    //     addPoint(x,y);
    //     coordinates.push_back(x);
    //     coordinates.push_back(y);
    // }

    auto x = (xmax+xmin)/2;
    auto y = (ymax+ymin)/2;
    addPoint(x, y);
    coordinates.push_back(x);
    coordinates.push_back(y);

    return 0; 
    // cout << "filled coordinates of (" << i()<<", "<<j()<<
    // ") on processor "<<r()<<endl;

}


void Bucket::sendCoordinates(int destination, int no_bucket) {

    // cout << r() <<": in send coordinates" << endl;
    auto ii = no_bucket/N;
    auto jj = no_bucket%N; 

    if (i() != ii || j() != jj ) {
        (*root)(ii,jj)->sendCoordinates(destination, no_bucket);
        return;
    }
    /* gerenate data if not available */
    if (coordinates.empty())
        fillCoordinates(0); 
    // cout << "cordinates to send: ";
    for (auto coord : coordinates) {
        cout << coord << " ";
    }
    cout << endl;
    int tmp = -1;
    MPI_Send(&tmp,1,MPI_INT,destination,TAG_NOTIFY,MPI_COMM_WORLD);//, &request);
    MPI_Send(coordinates.data(), MAX_PNTS, MPI_DOUBLE, destination, 
                TAG_DATA, MPI_COMM_WORLD);
}

int Bucket::getPoints(vector<Point>& pnts, int stat) {
    /* get or gerenate data if not available, 
    check if coordinates are empty because points can
    remain empty if no point in bucket */
    int ret_val = 0;
    if (coordinates.empty()) {
        ret_val = fillCoordinates(stat);
        if (!ret_val) pnts = points;
    } else pnts = points;
    return ret_val;
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
    bool open_request = false;
    for (auto i=0; i<size; i++)  
        MPI_Irecv( buf+i, 1, MPI_INT, i, 
                TAG_NOTIFY, MPI_COMM_WORLD, &requests[i] ); 
    bool startup = true;
    current_point_index = 0;
    int step = 0;
    if (startup) {
        // only on root bucket for now!
        step = initialize(step);
        if (step==0) startup = false;
        // test();
        // calculateDelauney();
    }
    while(running) { 
        
        MPI_Waitsome( size, requests, &ndone, indices, statuses ); 
        for (auto i=0; i<ndone; i++) { 
            auto j = indices[i]; 
            // printf( "Msg from %d with tag %d and value %d\n",  
            //         statuses[i].MPI_SOURCE,  
            //         statuses[i].MPI_TAG,
            //         buf[j]); 
            if (statuses[i].MPI_SOURCE != r() && buf[j]!=-1 ) {
                sendCoordinates(statuses[i].MPI_SOURCE, buf[j]);
            }
            // else { /* stop if receive msg from oneself */
            //     cout << "stop running " << endl;
            //     running = false;
            // }
            if (startup) {
                // only on root bucket for now!
                step = initialize(step);
                if (step==0) startup = false;
                // test();
                // calculateDelauney();
            }
            MPI_Irecv( buf+j, 1, MPI_INT, j, 
                    TAG_NOTIFY, MPI_COMM_WORLD, &requests[j]); 
        } 
    } 
}

/* return values:
   -------------- 
   0: completed
    
*/
int Bucket::initialize(int step) {

    vector<Point> pnts;

    /* always called from root processor, so there should be nothing to send */
    if (getPoints(pnts, 0)) throw runtime_error("initialize called from wrong processor");

    if (pnts.empty()) throw runtime_error("test point list empty");

    /* current (i,j) */
    auto ii = i();
    auto jj = j();

    if (step == 0)
        poly = move(Polygon(pnts[current_point_index]));

    /* build polygon for first time */ 
    vector<Point> tmp;
    auto add = [this](const Point& p) {poly.addPoint(p);};
    
    if (step == 1)
        goto resume;

    if (step > 1) 
        throw runtime_error("invalid step in initialize");

    for (init_dir_i=-1; init_dir_i<2; init_dir_i+=2 ) {
        for (init_dir_j=-1; init_dir_j<2; init_dir_j+=2 ) {
            init_incr = 1;
            while (true){                
                // if (!(i+incr*dir_i==N || i+incr*dir_i==0 || j+incr*dir_j==N || i+incr*dir_j==0) ){
                //     if (dir_i==-1 && dir_j==-1) A[i+incr*dir_i][j+incr*dir_j].pop_back();
                // } // just for testing!

                /* add ghost point if out of bounds */
                if (ii+init_incr*init_dir_i==N || ii+init_incr*init_dir_i==0 || jj+init_incr*init_dir_j==N || jj+init_incr*init_dir_j==0) {
                    // stop = 1;
                    add(Point(max(init_dir_i, 0), max(init_dir_j,0)));
                    break;
                }
                
                /* just for testing */
                // if ((*self)(5,5)->getPoints(tmp, 0)) {
                //     return 1;
                //     resume:
                //     (*self)(5,5)->getPoints(tmp, 1);
                // }

                if ((*self)(ii+init_incr*init_dir_i,jj+init_incr*init_dir_j)->getPoints(tmp, 0)) {
                    return 1;
                    resume:
                    (*self)(ii+init_incr*init_dir_i,jj+init_incr*init_dir_j)->getPoints(tmp, 1);
                }
                
                /* add points (see sketch) and go in diagonal direction if cell is empty*/
                if (!tmp.empty()) {
                    for_each(tmp.begin(), tmp.end(), add);
                    tmp.clear();
                    break;
                }
                else init_incr++;
            }
        }
    }
    poly.calculateVoronoi();

    cout << r() <<": center: " << poly.c << endl;
    cout << r() <<": poly: " << poly << endl;
    poly.printPoints(to_string(r()));
    cout << "****** Proc " << r() << " finished ******" << endl;

    return 0;
}


void Bucket::calculateDelauney(){ // test on one processor
    // auto it = 0; // iterator
    // while (it <poly.points.size()) {
    //     begin:
    //     cout << endl << endl << "it: " << it<<
    //         " voroni.size(): " << poly.voronoi.size() <<
    //         " points.size(): " << poly.points.size() << endl;
    //     auto vor_size = poly.voronoi.size();
    //     auto rad = poly.radii[it];
    //     auto pnt = poly.voronoi[it];
    //     /* distance of buckets that can contain candidates */
    //     int n = rad*N+1;
    //     /* indices of bucket containing pnt */
    //     int pi = pnt.x*N;
    //     int pj = pnt.y*N; 

    //     poly.V.clear();
    //     /* calculate Delauney neighbour candidates for pnt */
    //     for (auto di=-n; di<n+1; di++) 
    //         for (auto dj=-n; dj<n+1; dj++) 
    //             for (auto p : (*self)(pi+di, pj+dj)->getPoints()) 
    //                 /* add if in right range and not already in points or center*/
    //                 if (Point::dist(p,pnt) < rad 
    //                   && std::find(poly.points.begin(), poly.points.end(), p) == poly.points.end()
    //                   && p!=poly.c)
    //                     poly.V.push_back(p);

    //     if (poly.V.empty()) {
    //         cout << "V empty for " << pnt << endl;
    //         it++;
    //     } else {
    //         auto v = poly.V.front();
    //         /*
    //          TODO: check with paper again
    //          not sure about that part differs a bit from paper */ 
    //         int first = -1, last = -1;
    //         for (auto k=0; k<vor_size; k++) {
    //             /* search first and last index in half plane 
    //                 dist(c,p) <= dist(v,p) -> p in half plane
    //             */
    //             if (Point::dist(poly.voronoi[k],poly.c) <= Point::dist(poly.voronoi[k],v) &&
    //                 Point::dist(poly.voronoi[(k+1)%vor_size],poly.c) > Point::dist(poly.voronoi[(k+1)%vor_size],v)) 
    //                     last = k;
    //             if (Point::dist(poly.voronoi[k],poly.c) > Point::dist(poly.voronoi[k],v) &&
    //                 Point::dist(poly.voronoi[(k+1)%vor_size],poly.c) <= Point::dist(poly.voronoi[(k+1)%vor_size],v)) 
    //                     first = (k+1)%vor_size; // it is actually (first-1)%poly.size() !
    //         }
    //         cout <<  "first: " << first << " last: " << last << " " << endl;
    //         cout << "v: " << v << endl;
    //         // cout << " size voronoi: " << vor_size << 
    //         //     " size V: " << poly.V.size() << " for " << pnt << endl;

    //         if (first==-1 || last == -1) {
    //             it++;
    //             continue;
    //         }
    //         Point o1;
    //         /* ATTENTION: in paper it is last +1 and first-1*/ 
    //         if (!circumcenter(o1, static_cast<Point>(poly.points[(last+1)%poly.points.size()]), v, static_cast<Point>(poly.c))){
    //             if (Point::dist(v,poly.c) < Point::dist(poly.points[(last+1)%poly.points.size()], poly.c)) {
    //                 poly.points.erase(poly.points.begin()+(last+1)%poly.points.size());
    //                 poly.addPoint(v);
    //                 poly.calculateVoronoi();
    //                 cout << "here 1" << endl;
    //                 goto begin;

    //             } else {
    //                 it++;
    //                 cout << "here 2" << endl;
    //                 goto begin;
    //             }
    //         }
    //         Point o2; 
    //         if (!circumcenter(o2, static_cast<Point>(poly.points[first]), v, static_cast<Point>(poly.c))) {
    //             if (Point::dist(v,poly.c) < Point::dist(poly.points[first], poly.c)) {
    //                 poly.points.erase(poly.points.begin()+first);
    //                 poly.addPoint(v);
    //                 poly.calculateVoronoi();
    //                 cout << "here 3" << endl;

    //                 goto begin;
    //             } else {
    //                 it ++;
    //                 cout << "here 4" << endl;
    //                 goto begin;
    //             }
    //         }
            
    //         cout << endl << "poly before " << poly << endl;

    //         /* erase values */

    //         for (int k=(last+1)%vor_size; k!=first; k=(k+1)%vor_size) {
    //             cout << "vor" <<  k << endl;
    //             poly.voronoi.erase(poly.voronoi.begin()+k);
    //         }
    //         for (int k=(last+1)%vor_size; k!=first; k=(k+1)%vor_size)
    //             poly.radii.erase(poly.radii.begin()+k);
    //         auto tmp = (first-1)%vor_size; // max was not working
    //         if (tmp <0) tmp =0;
    //         for (int k=(last+1)%vor_size; k!=tmp; k=(k+1)%vor_size)
    //         {
    //             cout << "pnts" << k << endl;
    //             poly.points.erase(poly.points.begin()+k);
    //         }

    //         // poly.voronoi.erase(poly.voronoi.begin()+(last+1)%vor_size, poly.voronoi.begin()+first);
    //         // poly.radii.erase(poly.radii.begin()+(last+1)%vor_size, poly.radii.begin()+first);
    //         // poly.points.erase(poly.points.begin()+(last+1)%vor_size, poly.points.begin()+tmp);
    //         cout << endl << "poly after erase " << poly << endl;

    //         /* add new values */
    //         poly.addPoint(v);
    //         // cout << (last+1)%vor_size << endl;
    //         poly.voronoi.insert(poly.voronoi.begin()+(last)%vor_size, o2);
    //         poly.voronoi.insert(poly.voronoi.begin()+(last)%vor_size, o1);
    //         poly.radii.insert(poly.radii.begin()+(last)%vor_size,Point::dist(o1,poly.c));
    //         poly.radii.insert(poly.radii.begin()+(last)%vor_size,Point::dist(o2,poly.c));

    //         cout << endl << "poly after addition  " << poly << endl;

    //         // cout << "o1: " << o1 << " o2: " << o2 << endl;
    //         poly.printPoints(to_string(it));
    //         it++;
    //         // it=(last+1)%poly.points.size();
    //         // if (it == 0) break;
    //     }
                    
    // }
    // cout << "finished" << endl;
    // poly.printPoints(to_string(42));
    // printList();
    // /* calculate Delauney neighbour candidates  */
    // // poly.V.clear();
    // // for (auto it = 0; it < vor_size; ++it) {

    // //     auto rad = poly.radii[it];
    // //     auto pnt = poly.voronoi[it];
    // //     int n = rad*N+1;
    // //     int pi = pnt.x*N;
    // //     int pj = pnt.y*N; // indices of point 

    // //     for (auto di=-n; di<n+1; di++) 
    // //         for (auto dj=-n; dj<n+1; dj++) 
    // //             for (auto p : (*a)(pi+di, pj+dj)->getPoints()) 
    // //                 if (Point::dist(p,pnt) < rad) {
    // //                     poly.V.push_back(p);
    // //                     // poly.addPoint(p);
    // //                 }
    // // }

    // // poly.calculateVoronoi();
    // // cout << "poly: " << poly << endl;

    // // cout << "V: [ ";
    // // for (const auto& p: poly.V)
    // //     cout << p << " ";
    // // cout << "]" << endl;

}
/* test function to be deleted in the end */
void Bucket::test() { 

    // shared_ptr<Bucket> a = self;
    // a = (*a)(5,5);
    // // cout << "here " << endl;
    // // a = (*a)(min(i()+5, N-1),min(i()+5, N-1));  // this is a problem
    // // a = (*a)(N-1,N-1); 
    // // a = (*a)(5,5); 
    // // a = (*a)(5,5);
    // // a = (*a)(0,7);
    // // a = (*a)(1,N-2);
    // // a = (*a)(3,1);
    // // a = (*a)(7,7);
    // // problem: deadlock somewhere
    // for (int ii =5;ii<N;ii++) {
    //     for(int jj=5;jj<N;jj++) {
    //         a = (*a)(ii,jj);
    //         a->getPoints();
    //     }
    // }
    // // cout << "r: " << a->r() << endl;
    // vector<Point> pnts = move(a->getPoints());
    // // cout << "points on " << r() << " ";
    // // for (auto p: pnts)
    // //     cout << p << " ";
    // // cout << endl;

    // // int no_bucket = i()*N+j();
    // // MPI_Send(&no_bucket,1,MPI_INT,r(),TAG_NOTIFY,MPI_COMM_WORLD);


    // printList();
    // cout << "****** Proc " << r() << " finished ******" << endl;
    // // cout << "coordinates: ";
    // // for (auto p: coordinates)
    // //     cout << p << " ";
    // // cout << endl;

    // // a->printNeighbours();
    // // a->printList();

    // /* send stop notification */
    // // int stop = -1;
    // // MPI_Send(&stop,1,MPI_INT,r(),TAG_NOTIFY,MPI_COMM_WORLD);
    // cout << "here " << r() << endl;
}

/********* POLYGON *********/

template <typename T> 
bool circumcenter(T& ret, const T& a, const T& b, const T& c) {
    auto d = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if (d == 0)
        return false; 
        // throw runtime_error("lines are parallel in circumcenter");
    ret.x = ((a.x * a.x + a.y * a.y) * (b.y - c.y) + (b.x * b.x + b.y * b.y) * (c.y - a.y) + (c.x * c.x + c.y * c.y) * (a.y - b.y)) / d;
    ret.y = ((a.x * a.x + a.y * a.y) * (c.x - b.x) + (b.x * b.x + b.y * b.y) * (a.x - c.x) + (c.x * c.x + c.y * c.y) * (b.x - a.x)) / d;
    return true;
}
template bool circumcenter<Point>(Point& ret, const Point& a, const Point& b, const Point& c);

// add points in counter clockwise order
void Polygon::addPoint(double p0, double p1) {
    addPoint(PointPoly(this, p0, p1));
}

void Polygon::addPoint(const Point& pin) {
    PointPoly p(this, pin);
    // if (std::find(points.begin(), points.end(), p) != points.end())
    //     return; 
    bool exists = false;
    auto pos = std::find_if(points.begin(), points.end(), [p,&exists](auto s) {
        if (s==p) {
            exists=true;
        }
        return s < p;
    });
    if (exists) return; // if element alread exists don't add
    points.insert(pos, move(p));
}

void Polygon::calculateVoronoi(void) { // ToDo: spelling

    voronoi.clear();
    radii.clear();
    if (points.size() < 2) return;
    PointPoly cc(this); // to be circumcenter

    for(auto it = points.begin(); it!=points.end()-1; it++){
        auto next = it+1;
        circumcenter(cc, *it, *next, c);
        voronoi.push_back(cc); // cc is copied
        radii.push_back(Point::dist(cc,c));
    }
    circumcenter(cc, points.front(), points.back(), c);
    voronoi.push_back(cc);
    radii.push_back(Point::dist(cc,c));
}

std::ostream &operator<<(std::ostream &os, const Polygon &poly) {
    os << "points: [ ";
    for (const auto& p: poly.points)
        os << p << " ";
    os << "]" << endl << "voronoi: [ " ;
    for (const auto& p: poly.voronoi)
        os << p << " ";
    os << "]" << endl << "radii: [ " ;
    for (const auto& p: poly.radii)
        os << p << " ";
    os << "]";
    return os;
}

void Polygon::printPoints(string no) {

    ofstream myfile;
    
    myfile.open("polyPoints"+no+".txt");
    if (myfile.is_open()) {
        myfile << "[" << c << ",";

        for (const auto& p: points)
                myfile << p << ", ";
        myfile << "]" << endl;
        myfile.close();
    }
    else cout << "Unable to open file";

    myfile.open("voronoiPoints"+no+".txt");
    if (myfile.is_open()) {
        myfile << "[";
        for (const auto& p: voronoi)
                myfile << p << ", ";
        myfile << "]" << endl;
        myfile.close();
    }
    else cout << "Unable to open file";
}

bool Polygon::PointPoly::operator> (const Polygon::PointPoly& other) const {
    PointPoly a = *this-poly->c;
    PointPoly b = other-poly->c;

    // to have a beginning (at angle == 0 in unit circle)
    if (a.y>=0 && b.y<0) return 1;
    if (a.y<0 && b.y>=0) return 0;
    if (a.y==0 && b.y==0) {
        if (a.x>0 || b.x>0) return a.x<b.x;
        else return a.x>b.x;
    }

    // check order 
    auto cross_prod = a.x*b.y-a.y*b.x;
    if (cross_prod > 0) return 1;
    if (cross_prod < 0) return 0;

    // if points have same angle take nearer one first
    if (a.x*a.x+a.y*a.y < b.x*b.x+b.y*b.y) return 1;
    
    return 0;
}

bool Polygon::PointPoly::operator>= (const Polygon::PointPoly& other) const {
    if (*this>other || *this==other) return 1;
    return 0;
}

bool Polygon::PointPoly::operator< (const Polygon::PointPoly& other) const {
    if (!(*this>=other)) return 1;
    return 0;
}

bool Polygon::PointPoly::operator<= (const Polygon::PointPoly& other) const {
    if (!(*this>other)) return 1;
    return 0;
}
