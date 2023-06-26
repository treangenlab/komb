#ifndef KOMB_HASHINDEXEDMINHEAP_H
#define KOMB_HASHINDEXEDMINHEAP_H

#include <iostream>
#include <cstdlib>
#include <utility>
#include <climits>


class HashIndexedMinHeap
{
	public:
		HashIndexedMinHeap(int capacity)
		{
			this->array = (int*)malloc(sizeof(int)*capacity);
			this->positions = (int*)malloc(sizeof(int)*capacity);
			this->values = (double*)malloc(sizeof(double)*capacity);

			this->size = 0;
			this->capacity = capacity;

			for(int i = 0; i < this->capacity; i++)
			{
				this->positions[i] = this->missingPosition;
			}
		}

		~HashIndexedMinHeap()
		{
			free(array);
            		free(positions);
            		free(values);
		}

		int getSize()
		{
			return this->size;
		}

		bool containsKey(int key)
		{
			return (this->positions[key] == this->missingPosition) ? false : true;
		}

		std::pair<int,double> peek()
		{
        		if(size == 0)
        		{
            			return std::make_pair(INT_MIN,INT_MIN);
       			}
        
        		return std::make_pair(this->array[0], this->values[array[0]]);
    		}

		std::pair<int,double> poll()
		{
			if (this->size == 0)
			{
				return std::make_pair(INT_MIN,INT_MIN);
			}

			std::pair<int,double> top = peek();
			
			this->positions[top.first] = this->missingPosition;

        		if(size != 1)
        		{
            			int last = this->array[size-1];
            			this->array[0] = last;
            			this->positions[last] = 0;

            			this->size--;
            			this->minHeapfy(0);
        		}
        		else
        		{
            			this->size--;
        		}
        
        		this->array[size] = 0;

        		return top;
		}

    		bool insert(int key, double value)
    		{

        		if(this->size >= this->capacity)
        		{
            			return false;
        		}

        		int pos = this->size;
        		this->size++;
        		this->array[pos] = key;
        		this->positions[key] = pos;
        		this->values[key] = value;
        		this->refreshPriority(key, value);
        	
        		return true;
    		}

    		bool satisifiesHeap()
    		{

        		for(int pos = 0; pos < size; pos++)
        		{
            			int keyCur = this->array[pos];
            			int posLeft = (2*(pos+1))-1;
            			if(posLeft < this->size)
            			{
                			int keyLeft = this->array[posLeft];
                			if(this->values[keyCur] > this->values[keyLeft])
                			{
                    				return false;
                			}
            			}

            			int posRight = (2*(pos+1));
            			if(posRight < this->size)
            			{
                			int keyRight = this->array[posRight];
                			if(this->values[keyCur] > this->values[keyRight])
                			{
                    				return false;
                			}
            			}
        		}

        		return true;
    		}

    		double getPriority(int key)
    		{
 			return this->values[key];
    		}	

    		void refreshPriority(int key, double value)
    		{

        		this->values[key] = value;
        		int pos = this->positions[key];
        		bool shiftedDown = this->minHeapfy(pos);

        		if(!shiftedDown)
        		{
            			if(pos > 0)
            			{
        				int cur = key;
                			int parentPos = ((pos + 1) / 2) - 1;
                			int pel = this->array[parentPos];
                			while(pos > 0 && this->values[pel] > this->values[cur])
                			{
                    				this->array[parentPos] = cur;
                    				this->positions[cur] = parentPos;
                    				this->array[pos] = pel;
                    				this->positions[pel] = pos;
                    				pos = parentPos;
                    				parentPos = ((pos + 1) / 2) - 1;
                    				if(pos > 0)
                    				{
                        				pel = this->array[parentPos];
                    				}
                			}
            			}
        		}
    		}



		bool minHeapfy(int pos)
		{

        		int posLeft = (2*(pos+1))-1;
        		int posRight = (2*(pos+1));

        		int keyCur = this->array[pos];

        		int smallest = pos;
        		int nsmallest = keyCur;

        		if(posLeft < size)
        		{
            			int keyLeft = this->array[posLeft];
            			if(this->values[keyLeft] < this->values[keyCur])
            			{
                			smallest = posLeft;
                			nsmallest = keyLeft;
            			}
        		}

        		if(posRight < size)
        		{
            			int keyRight = this->array[posRight];
            			if(this->values[keyRight] < this->values[nsmallest])
            			{
                			smallest = posRight;
                			nsmallest = keyRight;
            			}
        		}

        		if(smallest != pos)
        		{

            			this->array[pos] = nsmallest;
            			this->positions[nsmallest] = pos;

            			this->array[smallest] = keyCur;
            			this->positions[keyCur] = smallest;

            			this->minHeapfy(smallest);
            			return true;
        		}

        		return false;

    		}




	private:

		/* array of keys */
		int* array;
		/* index -> position */
		int* positions;
		/*  index -> values */
		double* values;

		/* current size */
		int size;
		/* max capacity */
		int capacity;
		/* position to indicate missing key */
		const int missingPosition = -1;

};
#endif
