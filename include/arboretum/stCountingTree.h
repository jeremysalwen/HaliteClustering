// -*-c++-*-

#ifndef __STCOUNTINGTREE_H
#define __STCOUNTINGTREE_H

#include "stCell.h"
#include <db_cxx.h>
#include <string>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "Normalization.h"
namespace Halite {
  
  template<typename D>
  class stCountingTree {
  public:
    stCountingTree(int H, DBTYPE dbType, uint64_t cache_size, const std::string& tmpdir, int DIM) {
      //empty tree
      this->H = H;
      sumOfPoints = 0;

      P.resize(DIM,0);

      env=std::unique_ptr<DbEnv>(new DbEnv(0));

      u_int32_t cacheSizeG, cacheSizeB;
      cacheSizeG = cache_size / (1024*1024*1024);
      cacheSizeB = cache_size - (1024*1024*1024)*(uint64_t)cacheSizeG;
   
      if(env->set_cachesize(cacheSizeG, cacheSizeB, 1)) {
	std::cerr << "Error setting cache size.\n";
      }
      env->set_tmp_dir(tmpdir.c_str());
      
      if(env->open(NULL, DB_CREATE | DB_INIT_MPOOL | DB_PRIVATE, 0)) {
	std::cerr<< "Error opening db environment.\n";
      }
      for(int i=0; i<H-1; i++) {
	levels.push_back(std::unique_ptr<Db>(new Db(env.get(),0)));
      }
    
      //create/open one dataset per tree level
      for (int i=0; i<H-1; i++) {
	//create dbName
	std::string dbName=str(boost::format("level_%1%.db") % i);

	//boost::filesystem::remove(dbName);
	
	//create and open dataset
	levels[i]->set_flags(0); //no duplicates neither special configurations
	if(levels[i]->open(NULL, NULL, NULL, dbType, DB_CREATE, 0)) {
	  std::cerr<<"Error opening db.\n";
	}
      }
    }
    ~stCountingTree() {

      //close and remove one dataset per tree level
      for (size_t i=0; i<H-1; i++) {
	levels[i]->close(DB_NOSYNC); 
      }
      env->close(0);
      /*      for (size_t i=0; i<H-1; i++) {
  	std::string dbName=str(boost::format("level_%1%.db") % i);
	boost::filesystem::remove(dbName);
	}*/

    }

    char insertPoint(const D *point, D *normPoint) {
      const size_t nPos = (P.size() + 7) / 8;

      //creates used arrays
      std::vector<D> min(P.size(), 0);
      std::vector<D> max(P.size(), 1);
      std::vector<unsigned char> fullId((H-1)*nPos, 0);

      std::vector<stCell> cellAndParents;
      for (size_t i = 0; i < H - 1; ++i)
	cellAndParents.push_back(stCell(P.size()));

      // normalizes the point
      char out = 0;
      normalization->normalize(point, normPoint);
      for (unsigned int i=0; i<P.size(); i++) {
	if (normPoint[i] < 0 || normPoint[i] > 1) {
	  out = 1; // invalid point
	}
      }

      if (!out) { // if the point is valid
	sumOfPoints++; // counts this point in the root
	// inserts this point in all tree levels
	deepInsert(0, min.data(), max.data(), normPoint, fullId, cellAndParents);
      }//end if

      // return the point "validity"
      return !out;
    }

    Db *getDb(size_t level) {
      return (level<H-1)?levels[level].get():NULL;
    }
 
    int getSumOfPoints() {
      return sumOfPoints;
    }
    int *getP() {
      return &P[0];
    }

    void setNormalization(std::shared_ptr<Normalization<D>> normalization) {
      this->normalization=normalization;
    }
    int findInNode(std::vector<stCell>& parents, stCell* cell, stCellId& id, int level) {
      const size_t nPos = (P.size() + 7) / 8;
      std::vector<unsigned char> fullId((level+1)*nPos, 0);

      //use id and parents to define fullId
      for (int l=0; l<level; l++) {
	parents[l].getId()->getIndex(&fullId[l*nPos]);
      }
      id.getIndex(&fullId[level*nPos]);

      //search the cell in current level
      Dbt searchKey(fullId.data(),(level+1)*nPos), searchData;
      std::vector<unsigned char> cellSerialized (stCell::size(P.size()));
      searchData.set_data(cellSerialized.data());
      searchData.set_ulen(stCell::size(P.size())); //Dynamically finds the size of an stCell with dimension P.size()
      searchData.set_flags(DB_DBT_USERMEM);
      char found = (levels[level]->get(NULL, &searchKey, &searchData, 0) != DB_NOTFOUND);
      if (found) {
	*cell = stCell::deserialize(cellSerialized.data());
	cell->id = id;
      }

      return found;
    }

    void dbSync() {
      for (size_t i=0; i<H-1; i++) {
	levels[i]->sync(0); //flushes any cached information to disk
      }
    }
    void commitCell(std::vector<stCell>& parents, stCell *cell, size_t level) {
      const size_t nPos = (P.size() + 7) / 8;
      std::vector<unsigned char> fullId((level+1)*nPos, 0);

      //use cell and parents to define fullId
      for (size_t l=0; l<level; l++) {
	parents[l].getId()->getIndex(&fullId[l*nPos]);
      }
      cell->getId()->getIndex(&fullId[level*nPos]);
    
      // insert/update the dataset
      Dbt key(fullId.data(),(level+1)*nPos);
      unsigned char* const serialized = cell->serialize();
      Dbt data(serialized,stCell::size(P.size())); //Dynamically finds the size of an stCell with dimension P.size()
      levels[level]->put(NULL, &key, &data, 0);

      delete[] serialized;
    }
    void findParents(unsigned char *fullId, std::vector<stCell>& parentsVector, int level){
      const size_t nPos = (P.size() + 7) / 8;
      Dbt searchKey, searchData;
      searchKey.set_data(fullId);
      searchData.set_ulen(stCell::size(P.size())); //Dynamically finds the size of an stCell with dimension P.size()
      searchData.set_flags(DB_DBT_USERMEM);
      stCellId id(P.size());
      for (int l=0; l<level; l++) {
	//search the parent in level l
	std::vector<unsigned char> serialized(stCell::size(P.size()));
	searchKey.set_size((l+1)*nPos);
	searchData.set_data(serialized.data());
	levels[l]->get(NULL, &searchKey, &searchData, 0); //search
	parentsVector[l] =stCell::deserialize(serialized.data());
      
	id.setIndex(fullId+(l*nPos));
	parentsVector[l].setId(&id); //copy the id
      }
    }

  private:
    void deepInsert(size_t level, D* min, D *max, D *point, std::vector<unsigned char>& fullId, std::vector<stCell>& cellAndParents) {
      if (level < H-1) {
	// mounts the id of the cell that covers point in the current level
	// and stores in min/max this cell's lower and upper bounds
	D middle;
	stCellId cellId = stCellId(P.size());
	for (unsigned int i=0; i<P.size(); i++) {
	  middle = (max[i] - min[i]) / 2 + min[i];
	  if (point[i] > middle) { // if point[i] == middle, considers that the point
	    // belongs to the cell in the lower half regarding i
	    cellId.invertBit(i,P.size()); // turns the bit i to one
	    min[i] = middle;
	  } else {
	    // the bit i remains zero
	    max[i] = middle;
	  }
	}

	// updates the half space counts in level-1
	if (!level) { // updates in the root (level == 0)
	  for (unsigned int i=0; i<P.size(); i++) {
	    if (!cellId.getBitValue(i,P.size())) {
	      P[i]++; // counts the point, since it is in the root's lower half regarding i
	    }
	  }
	} else { // updates in the father (level > 0)
	  cellAndParents[level-1].insertPointPartial(&cellId, P.size());
	}

	// mounts the full id for the current cell
	int nPos = (P.size() - 1) / 8 + 1;
	cellId.getIndex(&fullId[level*nPos]);

	//searchs the cell in current level
	Dbt searchKey(fullId.data(),(level+1)*nPos), searchData;
	unsigned char* serialized = new unsigned char[stCell::size(P.size())];
	searchData.set_data(serialized);
	searchData.set_ulen(stCell::size(P.size()));//Dynamically finds the size of an stCell with dimension P.size()
	searchData.set_flags(DB_DBT_USERMEM);
	if (!levels[level]->get(NULL, &searchKey, &searchData, 0))
	  cellAndParents[level] = stCell::deserialize(serialized);
	delete[] serialized;

	//counts the new point
	cellAndParents[level].insertPoint();

	//inserts the point in the next level
	//and udates the partial counts for cellAndParents[level]
	deepInsert(level+1, min, max, point, fullId, cellAndParents);


	// insert/update the dataset
	Dbt key(fullId.data(),(level+1)*nPos);
	//stCell tmp= stCell(P.size());
	//cellAndParents[level].copy(&tmp,P.size());
	//tmp.reset();
	serialized = cellAndParents[level].serialize();
	Dbt data(serialized, stCell::size(P.size()));

	levels[level]->put(NULL, &key, &data, 0);
	delete[] serialized;
      }
    }

    std::vector<int> P;
    size_t sumOfPoints;
    size_t H;
    std::vector<std::unique_ptr<Db>> levels;
    std::unique_ptr<DbEnv> env;
    std::shared_ptr<Halite::Normalization<D>> normalization;
  };

}
#endif //__STCOUNTINGTREE_H
