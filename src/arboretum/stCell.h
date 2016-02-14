#ifndef __STCELL_H
#define __STCELL_H

#include "stCellId.h"

class stCell {
    public:
        static stCell* create(const size_t DIM)
        {
            return new stCell(DIM);
        }

        stCell() : usedCell(0), sumOfPoints(0), id(stCellId()) {
        }

        stCell(const size_t DIM) :
                usedCell(0), sumOfPoints(0), id(stCellId(DIM)) {
            P.resize(DIM);
        }
        
        stCell(const stCell& cell) :
                usedCell(cell.usedCell), sumOfPoints(cell.sumOfPoints), 
                id(cell.id), P(cell.P) {
        }
        
        stCell& operator= (const stCell &cell) {  
            usedCell = cell.usedCell;
            sumOfPoints = cell.sumOfPoints;
            id = cell.id;
            P = cell.P;
            return *this;
        }

        unsigned char* serialize() const {
            size_t size = this->size();
            unsigned char* result = new unsigned char[size];
            unsigned char* p = result;

            *(char*)p = usedCell;   p += sizeof(char);
            *(int*)p = sumOfPoints; p += sizeof(int);
            p += sizeof(stCellId);
            *(size_t*)p = P.size(); p += sizeof(size_t);

            //std::cout << "Serialized to size " << P.size() << std::endl << std::flush;
            for (size_t i = 0; i < P.size(); ++i) {
                *(size_t*)p = P[i]; p += sizeof(size_t);
            }

            return result;
        }

        size_t size() const {
            return sizeof(char) +
                sizeof(int) +
                sizeof(stCellId) +
                sizeof(size_t) * (P.size() + 1);
        }

        static size_t size(const size_t p) {
            return sizeof(char) +
                sizeof(int) +
                sizeof(stCellId) +
                sizeof(size_t) * (p + 1);
        }

        static stCell deserialize(unsigned char* s) {
            char usedCell = *(char*)s;  s += sizeof(char);
            int sumPoints = *(int*)s;   s += sizeof(int);
            s += sizeof(stCellId);
            size_t size = *(size_t*)s;  s += sizeof(size_t);
            std::vector<size_t> d;
            for (size_t i = 0; i < size; ++i) {
                d.push_back(*(size_t*)s); s += sizeof(size_t);
            }

            stCell result(size);
            result.usedCell = usedCell;
            result.sumOfPoints = sumPoints;
            result.P = d;

            return result;
        }

        void insertPoint() {
            sumOfPoints++;          
        }
        int getSumOfPoints() {
            return sumOfPoints;
        }
        char getUsedCell() {
            return usedCell;
        }   
        void useCell() {
            usedCell = 1;
        }
        int getP(const size_t i) const {
            return P[i];
        }
        stCellId *getId() {
            return &id;
        }
        void setId(const stCellId *id) {
            if(id) {
                this->id = *id;
            } else {
                this->id.reset();
            }
        }
        void insertPointPartial(stCellId *sonsCellId, const size_t DIM) {
            for (size_t i = 0; i < DIM; i++) {
                if (!sonsCellId->getBitValue(i, DIM)) {
                    P[i]++;
                }
            }
        }
        void copy(stCell *cell, const size_t DIM) const { 
            cell->usedCell = usedCell;
            cell->P = P;
            cell->sumOfPoints = sumOfPoints;
            cell->setId(&id, DIM);
        }

        void reset() {
            id.reset();
        }

    public:


        char usedCell;     
        int sumOfPoints;
        stCellId id;
        std::vector<size_t> P;
};

#endif //__STCELL_H

