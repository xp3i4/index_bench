#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>
#include <math.h>
#include <seqan/basic.h>
#include <bitset>

#define nl std::endl;

using namespace seqan;

const unsigned shapelength = 25;
const unsigned shapeweight = 18;
const unsigned blocklimit = 32;


typedef Iterator<String<Dna> >::Type TIter;
typedef Shape<Dna, MinimizerShape<shapelength> > TShape;
typedef Shape<Dna, UngappedShape<shapelength> > TShape_u;
typedef Shape<Dna, SimpleMShape> TMShape;

//typedef Index<StringSet<DnaString>, IndexQGram<MinimizerShape<shapelength, shapeweight>, OpenAddressing > > TIndex;
typedef Index<StringSet<DnaString>, IndexQGram<MinimizerShape<shapelength>, OpenAddressing > > TIndex;
typedef Index<StringSet<DnaString>, IndexQGram<UngappedShape<shapelength>, OpenAddressing > > TIndex_u;
typedef Index<StringSet<DnaString>, IndexQGram<SimpleMShape, OpenAddressing > > TMIndex;

typedef Value<TShape>::Type HValue;

int mTest1(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    TShape shape;
    TIndex index(reads);
    uint64_t sum = 0, p = 0;
    double time = sysTime();
    unsigned count = 0;
    std::cout << "mTest1(): " << std::endl;
    createQGramIndexDirOnly(index);
    std::cout << "    lenght Dir = " << length(index.dir) << " length Text = " << lengthSum(indexText(index)) << std::endl;
    std::cout << "    getDir start sysTime(): " << sysTime() - time << std::endl;
    //for (float m = 0.25; m <=1; m += 0.25)
    //{
    count = 0;
    for(uint64_t k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - shape.span + 1; j++)
        //for (uint64_t j = 0; j < 100; j++)
        {
            hashNext(shape, it + j);
            //std::cout << shape.hValue << " " << shape.XValue << std::endl;
            p = getDir(index, shape);
            sum += index.dir[p + 1] - index.dir[p];
            count++;
        }
    }
    sum=_getBodyCounth(sum);
    std::cout << "    sum = " << sum << " count " << count << std::endl;
    std::cout << "    getDir end sysTime(): " << sysTime() - time << std::endl;
    //}
    std::cout << "    End mTest1()" << std::endl;

    return 0;
}

int mTest2(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    TShape shape;
    TIndex index(reads);

    uint64_t sum = 0;
    double time = sysTime();
    std::cout << "mTest2(): " << std::endl;
    for (unsigned j = 24; j < 31; j++)
    for (unsigned k = 22; k < 23; k++)
    {
        std::cout << j << " " << k << std::endl;
        //resize(shape, 30, k);
        //resize(index.shape, 30, k);
        shape.span=j;
        shape.weight=k;
        index.shape.weight=k;
        index.shape.span=j;
        createQGramIndexDirOnly(index);
        std::cout << "    lenght Dir = " << length(index.dir) << " length Text = " << lengthSum(indexText(index)) << std::endl;
        std::cout << "    getDir start sysTime(): " << sysTime() - time << std::endl;
        for(uint64_t k = 0; k < length(genome); k++)
        {
            TIter it = begin(genome[k]);
            hashInit(shape, it);
            for (uint64_t j = 0; j < length(genome[k]) - shape.span + 1; j++)
            {
                hashNext(shape, it + j);
                sum = index.dir[getDir(index, shape)];
            }
        }
        sum=_getBodyCounth(sum);
        std::cout << "    sum = " << sum << std::endl;
        std::cout << "    getDir end sysTime(): " << sysTime() - time << std::endl;
        std::cout << "    End mTest1()" << std::endl;
    }

    return 0;
}

int mTest3(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    TShape shape, shape1;
    TIndex index(reads);
    uint64_t sum = 0, p = 0;
    double time = sysTime();
    std::cout << "mTest3(): " << std::endl;
    createQGramIndexDirOnly(index);
    std::cout << "    getDir start sysTime(): " << sysTime() - time << std::endl;
    std::cout << "    length Dir = " << length(index.dir) << std::endl;
    std::cout << "    length Text = " << lengthSum(indexText(index)) << std::endl;
    std::cout << "    length SA = " << length(index.sa) << std::endl;
    for(uint64_t k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
            p = getDir(index, shape);
            for (uint64_t n = _getBodyCounth(index.dir[p]); n < _getBodyCounth(index.dir[p + 1]); n++)
            {
                sum ^= _getSA_i2(index.sa[n]);
            }
        }
    }
    std::cout << "    sum = " << sum << std::endl;
    std::cout << "    getDir end sysTime(): " << sysTime() - time << std::endl;
    std::cout << "    End mTest1()" << std::endl;

    return 0;
}



int uTest(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    TShape_u t_shape;
    TIndex_u index(reads);
    unsigned kmerLength = t_shape.span;
    uint64_t sum=0, count=0, p = 0;
    double time = sysTime();
    std::cout << "uTest():\n";
    std::cout << "    fullDirLength " << _fullDirLength(index) << std::endl;

    indexCreate(index, FibreDir());
    std::cout << "    getBucket start sysTime(): " << sysTime() - time << std::endl;
    //for (float m = 0.25; m <= 1; m += 0.25)
    //{
    count=0;
    for(unsigned k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(t_shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - kmerLength + 1; j++)
        //for (uint64_t j = 0; j < 100; j++)
        {
            hashNext(t_shape, it + j);
            p = getBucket(index.bucketMap, t_shape.hValue);
            sum += index.dir[p + 1] - index.dir[p];
            count++;
        }
    }
    std::cout << "    sum = " << sum << " count = " << count << std::endl;
    std::cout << "    getBucket end sysTime(): "<< sysTime() - time<< std::endl;
    //}
    std::cout << "    End uTest()" << std::endl;
    return 0;
}

int uTest3(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    TShape_u t_shape;
    TIndex_u index(reads);
    unsigned kmerLength = t_shape.span;
    uint64_t sum=0, count=0, p = 0;
    double time = sysTime();
    std::cout << "uTest3():\n";
    std::cout << "    fullDirLength " << _fullDirLength(index) << std::endl;

    indexCreate(index, FibreSADir());
    std::cout << "    getSA start sysTime(): " << sysTime() - time<< std::endl;
    for(unsigned k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(t_shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - kmerLength + 1; j++)
        //for (uint64_t j = 0; j < 100; j++)
        {
            hashNext(t_shape, it + j);
            p = getBucket(index.bucketMap, t_shape.hValue);
            for (uint64_t k = index.dir[p]; k < index.dir[p+1]; k++)
                sum ^= index.sa[k].i2;
        }
    }
    std::cout << "    sum = " << sum << " count = " << count << std::endl;
    std::cout << "    getSA end sysTime(): "<< sysTime() - time<< std::endl;
    std::cout << "    End uTest()" << std::endl;
    return 0;
}

int umTest(StringSet<DnaString> & reads,  StringSet<DnaString> & genome)
{

    std::cout << "umTest() " << std::endl;
    TShape_u t_shape;
    TIndex_u index(reads);
    unsigned kmerLength = t_shape.span;
    uint64_t s=0, sum=0;

    TShape shape;
    TIndex index1(reads);

    double time = sysTime();
    indexCreate(index, FibreDir());
    std::cout << "        create 1 time " << sysTime() - time << std::endl;
    createQGramIndexDirOnly(index1);

    std::cout << "        h2y(shape, BucketMap<uint64_t>::EMPTY) = " << h2y(shape, BucketMap<uint64_t>::EMPTY) << std::endl;
    uint64_t v1, v2, p1, p2;
    for(unsigned k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(t_shape, it);
        hashInit(shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - kmerLength + 1; j++)
        {
            hashNext(t_shape, it + j);
            hashNext(shape, it + j);
            p1 = getBucket(index.bucketMap, t_shape.hValue);
            p2 = getDir(index1, shape);
            if (index.dir[p1+1] - index.dir[p1] != _getBodyCounth(index1.dir[p2 + 1]) - _getBodyCounth(index1.dir[p2]))
                std::cout << "    counth unequal " <<  index.dir[p1+1] - index.dir[p1] << " " << _getBodyCounth(index1.dir[p2 + 1]) - _getBodyCounth(index1.dir[p2]) << " " << shape.hValue<< std::endl;
            if ((uint64_t)t_shape.hValue - (uint64_t)shape.hValue)
                std::cout << "    hValue unequal " << j << " " <<  t_shape.hValue << " " <<  shape.hValue << std::endl;
            v1=h2y(shape, index.bucketMap.qgramCode[getBucket(index.bucketMap, t_shape.hValue)]);
            v2=_getBodyValue(index1.dir[getDir(index1, shape)]);
            if (v1 != v2)
                if (v1 ^ h2y(shape, BucketMap<uint64_t>::EMPTY) || v2 ^ _bitEmpty)
                    std::cout << "    YValue unequal " << k << " " << j << " t_shape.hValue = " << t_shape.hValue << " getDir(index1, shape) = " << getDir(index1, shape) << " XValue = " << shape.XValue << " YValue = " << shape.YValue<< " v2 = " << v2 << " " << length(index1.dir) << std::endl;
        }
    }
    std::cout << "    s = " << s << std::endl;
    std::cout << "    End umTest() " << std::endl;
    return 0;
}

void opKmerAnalysis5(StringSet<DnaString> & genome, StringSet<DnaString> & reads)
{

    TShape shape;
    TIndex_u index(genome);
    uint64_t mask  = 131072 - 1, dn, count=0, sum = 0;
    uint64_t anchor[mask + 1] = {1<<63};
    String<uint64_t> result;
    resize(result, length(reads));
    double time = sysTime();
    //std::cout << length(result) << std::endl;
    std::cout << "opAnalysis5() \n";

    //createQGramIndexDirOnly(index);
    indexCreate(index, FibreSADir());
    //_initAnchor(anchor, mask);
    for (uint64_t j = 0; j < length(reads); j++)
    {
    //std::cout << "done " << std::endl;
        TIter it = begin(reads[j]);
        hashInit(shape, it);
        for (uint64_t h = 0; h < length(anchor); h++)
            anchor[h] = 0;
        uint64_t x = 1;
        sum=0;
        for (uint64_t k = 0; k < (length(reads[j]) - shape.span + 1); k++)
        {
            hashNext(shape, it + k);
            dn = getBucket(index.bucketMap, shape.hValue);
            uint64_t pre = ~0;
            if(index.dir[dn+1] - index.dir[dn] < 32)
            {
                sum += index.dir[dn+1] - index.dir[dn];
                for (uint64_t n = index.dir[dn]; n < index.dir[dn + 1]; n++)
                {
                    if (index.sa[n].i2 - pre > 1000)
                    {
                        anchor[x++] = ((index.sa[n].i2- k) << 20) + k;
                        pre = index.sa[n].i2;
                    }
                }
            }
        }
        anchor[0] = anchor[1];
        uint64_t max1 = 0, max2 = 0, k1 = 0, k2 =0, countk = 0, max=0, c_b=0,
        start=0, ak=anchor[0], cbb=0, cbs=0, mask_cb = (1<<20) - 1,
        sb=0;
        //for (uint64_t k = 0; k < length(anchor); k++)
        std::sort(begin(anchor), begin(anchor) + x);

        //std::cout << std::endl << j << " ";
        
        for (uint64_t k = 1; k <= x; k++)
        {
            if (anchor[k]-ak < (1000<<20))
            {
                cbb++;
                        }
            else
            {
                std::sort(begin(anchor)+sb, begin(anchor)+k, 
                [&mask_cb](uint64_t &a, uint64_t &b){return (a & mask_cb) < (b & mask_cb);});
                for (uint64_t m = sb+1; m < k; m++)
                {
                    if(((anchor[m]-anchor[m-1]) & mask_cb) > shapelength) 
                        c_b += shapelength;
                    else
                    {
                        c_b += (anchor[m] - anchor[m-1]) & mask_cb; 
                    }
                }
                if (c_b > max)
                {
                    max = c_b;
                    start = sb;
                }
                sb = k;
                ak = anchor[k];
                cbb = 1;
                c_b = shapelength;
                cbs=0;
            }

        }
        result[j] = (anchor[start] >> 20);
        std::cout << j << " 1  " << result[j] << std::endl;
        //std::cout << "1 2  " << result[j] << std::endl;
            //std::cout << " " << j << " " << length(reads[j]) << " " << (anchor[start] >> 20) << " " << (k1 << 10) << " " << max1 << " " << k2 << " " << max2 << " " << countk / (float)length(anchor) << " " << sum << " " << (float)sum / length(reads[j]) << " " << anchor[start] << " " << start << std::endl;
        //std:: cout << k << " " << (anchor[k] & _AnchorMask_i1_) << std::endl; //<< " " << std::bitset<64>(anchor[k]) << std::endl;
        sum ^= max1;
        if (max1 == 0)
            count++;
    }
    std::cout << "    End mnKmerAnalysis() systime() = " << sysTime() - time << std::endl;
}

int main(int argc, char** argv)
{
    if (argc < 3)
        return 1;
    //time = sysTime();
    SeqFileIn gFile(toCString(argv[1]));
    SeqFileIn rFile(toCString(argv[2]));
    StringSet<CharString> ids_r;
    StringSet<CharString> ids_g;
    StringSet<DnaString> reads;
    StringSet<DnaString> genome;
    readRecords(ids_r, reads, rFile);
    readRecords(ids_g, genome, gFile);
    //std::cout << "read done sysTime " << sysTime() - time << std::endl;
    //umTest(reads, genome);
    //uTest(reads, genome);
    //mTest1(reads, genome);
    //mTest2(reads, genome);
    //uTest3(reads, genome);
    //mTest3(reads, genome);
    //mTest4(reads, genome);
    opKmerAnalysis5(genome, reads);
    return 0;
}
