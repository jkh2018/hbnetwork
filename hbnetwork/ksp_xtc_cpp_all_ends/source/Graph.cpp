
#include "Graph.h"

using namespace Graph_all;

#define INF 9999999


Graph::Graph() : _N(0), _size(0)
{
}

Graph::Graph(dag_t& array)
{
	_N = _size = array.size();
	_array.resize(_size);
	for(unsigned int i=0; i<_size; ++i)
	{
		_array[i].resize(_size);
		for(unsigned int j=0; j<_size; ++j)
			_array[i][j] = array[i][j];
	}
}

Graph::Graph(double **array, unsigned N) : _N(N), _size(N)
{	
	_array.resize(_size);
	for(unsigned int i=0;i<_size;++i)
	{
		_array[i].resize(_size);
		for(unsigned int j=0;j<_size;++j)
			_array[i][j] = array[i][j];
	}
	
//	this->Output(); // for debug
}


Graph::~Graph()
{	
}

void Graph::Restart(dag_t& array)
{
	_array.clear();

	_N = _size = array.size();
	_array.resize(_size);
	for(unsigned int i=0;i<_size;++i)
	{
		_array[i].resize(_size);
		for(unsigned int j=0;j<_size;++j)
			_array[i][j] = array[i][j];
	}

//	this->Output(); // for debug
}

void Graph::Restart(double** array,unsigned N)
{
    	_array.clear();

	_N = _size = N;
	_array.resize(_size);
	for(unsigned int i=0;i<_size;++i)
	{
		_array[i].resize(_size);
		for(unsigned int j=0;j<_size;++j)
			_array[i][j] = array[i][j];
	}
}

double Graph::Dijkstra(unsigned start , unsigned end , path_t& path )
{
	if(end == start || start >= _size) return -1;

	int* paths = new int[_size];
	std::cout << "process OK 1" << std::endl;
	double min = dijkstra(paths,start,end);
	std::cout << "process OK 2" << std::endl;
        if(min < 0) { delete[] paths; return -1;}
	std::cout << "process OK 3" << std::endl;
	// parse the shortest path
	if(end > _size-1) end = _size - 1;
	int i = end;	
	std::cout << "process OK 4" << std::endl;
	while(i>=0)
	{		
		path.push(i);
		i=paths[i];			
	}
	// push the start node
// 	if(_array[start][path.top()]>=0)
// 		path.push(start);
// 	else 
// 		{ delete[] paths; return -1;}
//         std::cout << "last : " << i << " path top" << path.top() << std::endl;
	if(path.top() != start)
        {
            delete[] paths; return -1;
        }

        delete[] paths;
	return min;
}

void Graph::Output(ostream& out)
{
	for(unsigned int i=0;i<_size;++i)
	{
		out.width(2);
		out << std::right << i << " ";
		for(unsigned int j=0;j<_size;++j)
		{
			out.width(10); //out.precision(8);
			out << std::right <<_array[i][j] <<" ";
		}
		out << std::endl;
	}
}


/*
    S       -- START node.
	T		-- END node
			out.width(10); //out.precision(8);
			out << std::right <<_array[i][j] <<" ";
		}
		out << std::endl;
	}
}
*/
/*
    S       -- START node.
	T		-- END node
	Path[]  -- Path[i]=j indicates the previous node of node i is j. Path[i]=-1 indicates node i
	           has not previous node, e.g. START node. It could be expanded dynamically.
        Dist[]  -- store the lengths of the shortest pathes from START node to END node. It could be
			   expanded dynamically.Dist[i]=INF indicates there is not a path from START
			   node to node i.
	Prime[] -- Prime[i] indicates node Prime[i] is the prime node of node i. Prime[i]=-1 indicates
	           the node i has not prime node.// artificial node of i, which has broken connection with previous node.
	           To calculate kth path, we break nearest neighbors to save new node(_size+1, ...+2) and save Prime[i]=new
	Base[]  -- In opposition to Prime[i], node Base[i] is the original node of node i, that is to say, 
	           node i is the prime node of node Base[i]. //original node of prime node
	
 */

//KSP between start and end nodes, all data are similar to the previous one
int Graph::MSAforKSP(unsigned start, unsigned end ,unsigned k,  path_list& kpaths, int* path_length )
{
	if(end == start || start >= _size) return -1;
	// std::cout << " start " << start << std::endl;
	// std::cout << " end   " << end   << std::endl;
	// std::cout << " size " << _size << std::endl;
	/* Initialize */
	vector<int> Path, Prime, Base;
	vector<double> Dist;
	Path.resize(_size);
	Prime.resize(_size);
	Base.resize(_size);
	Dist.resize(_size);

	// std::cout << " end of initialization of k-path" << std::endl;

	for(unsigned int i=0; i<_size; ++i)
	{ Prime[i] = -1; Base[i] = i; }
	
	/* Find the shortest path firstly */
	int min = dijkstra(start, end, &Path[0], &Dist[0]);
	// std::cout << "shortest path (dijkstra) calculation finished!" << std::endl;
	// std::cout << "shortest value " << min << "INF " << INF << std::endl;
	if( min >=INF ) return 0;
	
	path_t path;
// 	std::cout << "shortest path value finished!" << std::endl;
	if(end > _size-1) end = _size - 1;
	int j=end;
	int ttt = 0;
	while(j>=0)
	{
		path.push(j); j=Path[j]; ttt += 1;
// 		std::cout << "path node index" << j << "path length " << ttt << std::endl;
	}
	path_length[0]=ttt;
	
// 	// push the start node
// 	if(_array[start][path.top()]>=0)
// 		path.push(start);
//         std::cout << "last : " << j << " path top" << path.top() << std::endl;
	if(path.top() != start)
        {
            return -1;
        }
	kpaths.push_back(path); // store the shortest path

//	 Find the 2th - kth shortest paths
	path.pop(); //path includes start node, but in real calculation,
				// we pick up the point next to start,
				// so we cut the start at the starting point.
	// int min_path_length = path_length[0];
	
	unsigned int ki = 1;
// 	while( (ki < k) && (path_length[ki-1] == min_path_length) )
	while(ki < k)
	{
		/* Find the first node with more than a single incoming arc */
		int nh = -1;
//                 std::cout << "path size " << path.size() << " ki " << ki << std::endl;
		while( path.size() )
		{
			unsigned node = path.top(); path.pop();
			int count = 0;
			for(unsigned int i=0; i<_size; ++i )
			{
				if( (_array[i][node] >= 0) && (i != start) ) count++;//node is next to start, so if
                                                                                    //other points to end into node, we pick up
				if( count == 1 ) break;
			}

			if( count == 1 ) { nh = node; break; }
		}

		if( nh < 0 ) break; // there is NOT an alternative path, exit!

		int ni = -1;
		/* Add the first prime node to graph */
		if( Prime[nh] < 0 )
		{
			unsigned nh1 = addNode(nh,Path[nh]); //nh1 is new(temparary) _size+1 node which has the same connection as nh
                                                            // but Path[nh], i.e., previous node of nh on the path is not connected with nh

			/* compute the minimal distance from node start to nh1 */
			double min_dist = (unsigned)INF;
			int min_node = -1;
			for(unsigned int i=0;i<_size;++i)
			{
// 				if(i == end) continue;
// 				std::cout << "i=" << i << "nh " << nh << "Prev: " << Path[nh] <<" : nh1 " << nh1 << " " <<Dist[i] << " " << _array[i][nh1] << std::endl; // for debug
				if( Dist[i] + (unsigned)_array[i][nh1] < min_dist )
				{
					min_dist = Dist[i] + _array[i][nh1];
					min_node = i;
				}
			}
//                         std::cout << "min_dist" << min_dist << "min_node" << min_node << std::endl;
			Dist.push_back(min_dist);
			Path.push_back(min_node);
// 			std::cout << min_dist << "min_node" << min_node << std::endl;
			Prime.push_back(-1); //Prime[nh1]=-1 because nh1 is artifical one
			Prime[nh] = nh1;

			/* record the base node */
			unsigned basei = nh;
			while(basei != (unsigned)Base[basei])
				basei = Base[basei];
			Base.push_back(basei); //Base(nh1)=nh
                        
//                         std::cout << "path size " << path.size() << " ki " << ki << std::endl;
                        
			if(path.size())
			{ ni = path.top(); path.pop(); }  //next node after nh on the path
		}
			/*  Get node ni, it must meet it's the first node following nh in path, but its prime node ni is NOT in graph */
		else
		{
			while( path.size() )
			{
				ni = path.top(); path.pop();
				if(Prime[ni] < 0) break;
			}
		}

		/* Add the other prime nodes to graph */
//                 std::cout << "ni " << ni <<std::endl;
		while(ni >=0 )
		{
			unsigned ni1 = addNode(ni,Path[ni]);
			if(Prime[Path[ni]]>=0)
				_array[Prime[Path[ni]]][ni1] = _array[Path[ni]][ni];	// add the arc -- (ni-1,ni)
									//For example, nh and ni are neighbors
									// nh and nh1 is pretty much the same but prev[nh] is not connected with nh1
									//For nh1, we set nh1 and ni is not connected, so we loose the connection from start to nh1
									// and then conntected into ni1, so we set nh1 and ni1 as the same as nh and ni

			/* compute the minimal distance from node start to ni1 */
			double min_dist = (unsigned)INF;
			int min_node = -1;
			for(unsigned int i=0;i<_size;++i)
			{
// 				if(i == end) continue;
// 				std::cout << "i=" << i << "ni " << ni << " Prev :" << Path[ni] <<" ni1: " << ni1 << " " << Dist[i] << " " << _array[i][ni1] << std::endl; // for debug
				if( Dist[i] + (unsigned)_array[i][ni1] < min_dist )
				{
					min_dist = Dist[i] + _array[i][ni1];
					min_node = i;
				}
			}

			Dist.push_back(min_dist);
			Path.push_back(min_node);
// 			std::cout << "min_dist" << min_dist << "min_node" << min_node << std::endl;
			Prime.push_back(-1);
			Prime[ni] = ni1;     //ni1 is prime of ni

			/* record the base node */
			unsigned basei = ni;
			while(basei != (unsigned) Base[basei])
				basei = Base[basei];
			Base.push_back(basei);  //set base of ni1 as ni, because of final result

			if( !path.size() ) break;
			ni = path.top(); path.pop();
		}

		/* get the kth shortest path */
		if( ni < 0 ) {ni = nh;}; // if nh is just the end node.
// 		std::cout << "what is ni" << ni << ", prime: " << Prime[ni] << ", base: " << Base[ni]<< std::endl;
		path_t temp; //kth shortest path file
		int j = Prime[ni];
		int ttt = 0;
		while(j>=0)
		{
			path.push(j); temp.push(Base[j]); j=Path[j]; ttt += 1;
		}
		path_length[ki] = ttt;
//                 std::cout << "temp size " << temp.size() << " ki " << ki << std::endl;
                if(temp.size()<2) break;

// 		// push the start node into the temp
// 		if(_array[start][temp.top()]>=0)
// 			temp.push(start);
//                 std::cout << "last : " << j << " path top" << path.top() << std::endl;
                if(temp.top() != start) break;

		kpaths.push_back(temp);  // store the kth shortest path

		ki++;
// 		this->Output(); // for debug
	}

	return ki;
}


/*
	Look out: Dijkstra algorithm doesn't assure the
	generated tree from graph is minimal generated tree.
*/
/// Look out: the returned paths doesn't include the start node!
double Graph::dijkstra( int* paths , unsigned start , unsigned end )
{	
	int* Used = new int[_size];		   // label the used nodes, 0 - indicates not used.
	double* Dist = new double[_size];     // store the min distances from start node to node indicated by current suffix.
//	std::cout << "real procedure of dijkstra" << std::endl;
	/* Initialize */	
	unsigned i;
	for(i=0;i<_size;++i)
	{
// 		paths[i]= -1;	
		paths[i]= (_array[start][i] == INF) ? -1 : start;
		Used[i] = 0;
		Dist[i] = _array[start][i]; 
	}
//	std::cout << "real procedure of dijkstra after intial" << std::endl;
//	std::cout << " start " << start << std::endl;
//	std::cout << " end   " << end << std::endl;
		
	Used[start] = 1;		           // label the start node
// 	Dist[start] = 0;
	
	unsigned count = 0;	
	while(count++ < _size)
	{		
		int min_node = -1;
		double min_dist = (unsigned)INF;

		/* Select */
		for(i=0; i<_size; ++i)
		{
			if(Used[i] > 0)continue;   // ignore the used node
			
			if(Dist[i] < min_dist)
			{
				min_dist = Dist[i];
				min_node = i;
//				std::cout << "min_node " << i << std::endl;
			}
		}
		if(min_dist == INF)break;
//		std::cout << "until ok!" << std::endl;
//		std::cout << "infinity " << INF << std::endl;
//		std::cout << " min_dist " << min_dist << std::endl;
		/* Record shortest path from the current node to start node */
		if(min_node >= 0)
		{			
			Used[min_node]  = 1;
			Dist[min_node]  = min_dist;			
		}
		if(min_node >= (int)_size-1)break;
//		std::cout << "until ok _1 !" << min_node << std::endl;
		/* Adjust */
		for(i=0; i<_size; ++i)
		{
			if( Used[i] > 0)continue;  // ignore the used node

			double w = (unsigned)_array[min_node][i]; //_array[min_node][i] < 0.0 ? INF : _array[min_node][i];
			if( min_dist + w < Dist[i] )
			{
				Dist[i] = min_dist + w;
				paths[i] = min_node;
			}
		}       
	}	
//	std::cout << "until ok _2 ! " << std::endl;
	double mdist = INF;
        mdist = Dist[_size-1];
	if( end != start && end < _N)
		mdist = Dist[end];
 //       std::cout << "mindist" << mdist << std::endl;
	delete[] Used;
	delete[] Dist;

	return (signed) mdist; // == INF ? -1 : mdist;
}

/*
	Look out: Dijkstra algorithm doesn't assure the
	generated tree from graph is minimal generated tree.
*/
/// Look out: the returned paths doesn't include the start node!
double Graph::dijkstra(unsigned start, unsigned end, int* paths, double* dists )
{
	int* Used = new int[_size];		   // label the used nodes, 0 - indicates not used.
	double* Dist = dists;     // store the min distances from start node to node indicated by current suffix.
                
	/* Initialize */
	unsigned int i;
	for(i=0;i<_size;++i)
	{
		Used[i] = 0;
//                 paths[i] = -1;
		paths[i]= (_array[start][i] >=INF ) ? -1 : start;
		Dist[i] = _array[start][i];
                
	}
	Used[start] = 1;		           // label the start node
        
	unsigned count = 0;
	while((count++ <= _size) && (Used[end] == 0))
	{
		int min_node = -1;
		double min_dist = (unsigned)INF;

		/* Select */
		for(i=0;i<_size;++i)
		{
			if(Used[i] > 0)continue;   // ignore the used node

			if(Dist[i] < min_dist)
			{
				min_dist = Dist[i];
				min_node = i;
			}
		}
  		//std::cout << "min_node " << min_node << " min_dist " << min_dist << "used end " << Used[end] << " count" << count << std::endl;
		if(min_dist >= INF)break;
                                
		/* Record shortest path from the current node to start node */
		if(min_node >= 0)
		{
			Used[min_node]  = 1;
			Dist[min_node]  = min_dist;
		}
		if(min_node >= (int)_size)break;
                //std::cout << "min_node " << min_node << " min_dist " << min_dist << "used end " << Used[end] << " count" << count << std::endl;
                
		/* Adjust */
		for(i=0; i<_size; ++i)
		{
			if( Used[i] > 0)continue;  // ignore the used node

			double w = (unsigned)_array[min_node][i]; //_array[min_node][i] < 0.0 ? INF : _array[min_node][i];
			if( min_dist + w < Dist[i] )
			{
				Dist[i] = min_dist + w;
				paths[i] = min_node;
			}
		}
	}

	//std::cout << Dist[end];
	double mdist = INF;
	if( end != start && end < _size)
		mdist = Dist[end];

	delete[] Used;

	return mdist; //== INF ? -1 : mdist;
}

//KSP between start and end nodes, all data are similar to the previous one
int Graph::dijkstra_all_ends(int start, int* end ,int num_ends,  path_list& kpaths, int* path_length )
{
        /* Initialize */
	unsigned int i;
	if(start >= _size) return -1;
        for(i=0; i<num_ends; ++i) { if(end[i] >= _size) return -1; }
	
	/* Find the shortest path firstly */
	int* Used = new int[_size];		   // label the used nodes, 0 - indicates not used.
	int* paths = new int[_size];
        int* Dist = new int[_size];
                
	
	for(i=0;i<_size;++i)
	{
		Used[i] = 0;
//                 paths[i] = -1;
		paths[i]= (_array[start][i] >=INF ) ? -1 : start;
		Dist[i] = _array[start][i];
                
	}
	Used[start] = 1;		           // label the start node
        
	unsigned count = 0;
        unsigned end_count = 0;

        while((count++ <= _size) && ( end_count <= num_ends))
	{
		int min_node = -1;
		double min_dist = (unsigned)INF;

		/* Select */
		for(i=0;i<_size;++i)
		{
			if(Used[i] > 0)continue;   // ignore the used node

			if(Dist[i] < min_dist)
			{
				min_dist = Dist[i];
				min_node = i;
			}
		}
  		//std::cout << "min_node " << min_node << " min_dist " << min_dist << "used end " << Used[end] << " count" << count << std::endl;
		if(min_dist >= INF)break;
                                
		/* Record shortest path from the current node to start node */
		if(min_node >= 0)
		{
			Used[min_node]  = 1;
			Dist[min_node]  = min_dist;
		}
		if(min_node >= (int)_size)break;
                //std::cout << "min_node " << min_node << " min_dist " << min_dist << "used end " << Used[end] << " count" << count << std::endl;
                
		/* Adjust */
		for(i=0; i<_size; ++i)
		{
			if( Used[i] > 0)continue;  // ignore the used node

			double w = (unsigned)_array[min_node][i]; //_array[min_node][i] < 0.0 ? INF : _array[min_node][i];
			if( min_dist + w < Dist[i] )
			{
				Dist[i] = min_dist + w;
				paths[i] = min_node;
			}
		}
		for(i=0; i<num_ends; ++i) { if(end[i] == min_node) end_count += 1; }
	}

	
// 	std::cout << "shortest path value finished!" << std::endl;
        for(i=0; i<num_ends; ++i) {
            path_t path;
            if(end[i] > _size-1) end[i] = _size - 1;
            int j=end[i];
            int ttt = 0;
            while(j>=0)
            {
		path.push(j); j=paths[j]; ttt += 1;
// 		std::cout << "path node index" << j << "path length " << ttt << std::endl;
            }
            path_length[i]=ttt;
	
            //if(path.top() != start)
            //{
            //    return -1;
            //}
            kpaths.push_back(path); // store the shortest path
        }

	return 0;
}

/*
	Here, using vector is more efficient than 2-dimension dynamic array.
*/
unsigned Graph::addNode(unsigned ni,int preni) //preni is the previous node of ni
{
	vector<double> newRow;
	for(unsigned int i=0;i<_size;++i)
	{
		if(i != (unsigned) preni)
			_array[i].push_back(_array[i][ni]);//Add new element [i][ni] into array[i] coloumn
		else
			_array[i].push_back(-1); //Delete connection between [preni][ni] for preni element of array[i] coloumn

		newRow.push_back(-1); //new vector _size elements of _size+1 sqaure matrix
	}
	newRow.push_back(-1); //Add one more element
	_array.push_back(newRow); //Push (Add) new row line

	return _size++; //size is one increased
}

vector<int> Graph :: prims_MST()
{
    vector<int> selected(_size); //selected nodes to be explored
    vector<int> tree_selected(_size); //selected nodes of individual tree in the case of several tree clusters
    vector<int> ndx_tree_nodes(_size); //tree index of nodes
    int n_search_node; //number of searched nodes
    int ndx_tree_const=0;     //index of tree
    int i,j;
    int min,x,y;
    
    //std::cout << " # of nodes " << nodes << endl;
    for(unsigned int i=0; i<_size; i++)
    {
        selected[i]=false;
        tree_selected[i] = false;
    }
    selected[0] = true;
    ndx_tree_nodes[0] = ndx_tree_const;
    n_search_node = 1;
    while( n_search_node <= _size-1 )
    {
        min=INF;
        for(unsigned int i=0; i<_size; ++i)
        {
            if((selected[i] == true) && (tree_selected[i] == false)){
                for(unsigned int j=0; j<_size; ++j){
                    if((selected[j] == false) && (tree_selected[i] == false)){
                        if(min > _array[i][j])
                        {
                            min=_array[i][j];
                            x=i;
                            y=j;
                        }
                    }
                }
            }
        }
        if(min == INF)
        {
            
            ndx_tree_const += 1;
            for(unsigned int i=0; i<_size; ++i)
            {
                if(selected[i] == true)
                {
                    tree_selected[i] = true;
                }
            }
            for(unsigned int i=0; i<_size; ++i)
            {
                if(selected[i] == false)
                { 
                    y = i;
                    //std::cout << " min values is INF " << " y " << y << " ndx_tree " << ndx_tree_const << endl;
                    break;
                }
            }
        }
        selected[y] = true;
        ndx_tree_nodes[y] = ndx_tree_const;
        n_search_node += 1;
        //std::cout <<x<<" --> "<<y << " " << ndx_tree_nodes[y] << endl;
    }
    
    return ndx_tree_nodes;
}
