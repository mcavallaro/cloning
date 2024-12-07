template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}

template<class T>
inline const T &MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
        {return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
        {return b > a ? float(b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
        {return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
        {return b < a ? float(b) : (a);}



//se voglio tenere in memoria i numeri piu' grandi, allora devo usare add_large e hs.larger(i)
//al contrario, se mi interessano i numeri piu piccoli allra devo usare add_small e hs.smaller(i)

template<class T>
void sort(std::vector<T> &arr,  std::vector<int> &index,  int m=-1){
	static const int M=7, NSTACK=64;
	int i,ir,j,k,jstack=-1,l=0,n=arr.size();
	if ( n != index.size() ){
		exit(1);
	}
	T a;
	int ai;
	std::vector<int> istack(NSTACK);
	if (m>0) n = MIN(m,n);
	ir=n-1;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				ai=index[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					index[i+1]=index[i];
				}
				arr[i+1]=a;
				index[i+1]=ai;
			}
			if (jstack < 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		}
		else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1]);
			SWAP(index[k],index[l+1]);
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir]);
				SWAP(index[l],index[ir]);
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir]);
				SWAP(index[l+1],index[ir]);
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1]);
				SWAP(index[l],index[l+1]);
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			ai=index[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
				SWAP(index[i],index[j]);
			}
			arr[l+1]=arr[j];
			index[l+1]=index[j];
			arr[j]=a;
			index[j]=ai;
			jstack += 2;
			if (jstack >= NSTACK)
				throw("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}



class Heapselect {
    private:
	int m,n,srtd;

    public:
	std::vector<double> heap;
	std::vector<int> indexes;
	Heapselect(int mm) : m(mm), n(0), srtd(0), heap(mm,1.e99), indexes(mm,mm+1) {
//        std::cout << "Heap is constructed" << std::endl;
    }

	~Heapselect() {
 //       std::cout << "Heap is destructed" << std::endl;
    }

	void add_large(double val, int index) {
		int j,k;
		if (n<m) {
			heap[n] = val;
			indexes[n] = index;
			n++;
			if (n==m) sort(heap,indexes);
		}
		else {
			if (val > heap[0]) {
				heap[0]=val;
				indexes[0] = index;
				for (j=0;;) {
					k=(j << 1) + 1;
					if (k > m-1)
						break;
					if (k != (m-1) && heap[k] > heap[k+1])
						k++;
					if (heap[j] <= heap[k])
						break;
					SWAP(heap[k],heap[j]);
					SWAP(indexes[k],indexes[j]);
					j=k;
				}
			}
			n++;
		}
		srtd = 0;
	}


	void add_small(double val, int index) {
		int j,k;
		if (n<m) {
			heap[n] = val;
			indexes[n] = index;
			n++;
			if (n==m) sort(heap,indexes);
		}
		else {
			if (val < heap[m-1]) {
				heap[m-1]=val;
				indexes[m-1] = index;
				for (j=0;;) {
					k=(j << 1) + 1;
					if (k > m-1)
						break;
					if (k != (m-1) && heap[k] > heap[k+1])
						k++;
					if (heap[j] <= heap[k])
						break;
					SWAP(heap[k],heap[j]);
					SWAP(indexes[k],indexes[j]);
					j=k;
				}
			}
			n++;
		}
		srtd = 0;
	}


	double smaller_value(int k){ //return  the k-th smallest value
		int mm = MIN(n,m);
		if (k >= mm) throw("Heapselect k too big");
		if (! srtd) {
			sort(heap,indexes);
			srtd = 1;
		}
		return heap[k];
	}


	double smaller_index(int k){ //return  the k-th smallest value
		int mm = MIN(n,m);
		if (k >= mm) throw("Heapselect k too big");
		if (! srtd) {
			sort(heap,indexes);
			srtd = 1;
		}
		return indexes[k];
	}


	double larger_index(int k) {//return  the k-th largest value	
		int mm = MIN(n,m);
		if (k >= mm) throw("Heapselect k too big");
		if (k == m-1) return heap[0];
		if (! srtd) {
			sort(heap,indexes);
			srtd = 1;
		}
		return indexes[mm-1-k];
	}

	double larger_value(int k) {//return  the k-th largest value	
		int mm = MIN(n,m);
		if (k >= mm) throw("Heapselect k too big");
		if (k == m-1) return heap[0];
		if (! srtd) {
			sort(heap,indexes);
			srtd = 1;
		}
		return heap[mm-1-k];
	}

};
