#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <map>
#include <vector>
#include <algorithm> 


using namespace std;
using namespace NTL;



void PopulateRandomMatrix(Mat<ZZ_p> &gfMatrix){
    for (size_t i = 0; i < gfMatrix.NumRows(); i++)
    {
        for (size_t j = 0; j < gfMatrix.NumCols(); j++)
        {
            ZZ_p randomNumber;
            random(randomNumber);
            gfMatrix.put(i, j, randomNumber);
        }

    }
}

ZZ_p GetRandomVWithP(Vec<float> distrib){
    int mod = conv<int>(ZZ_p::modulus());

    auto r = ((double) rand() / (RAND_MAX));
    int z = 0;
    for (z; z < mod; z++){
        if(r <= distrib[z])
            break;

    }
    return ZZ_p{z};
}


std::map<ZZ, float>  GetFrequency(Vec<ZZ_p> EntryVector){
    float entropy=0;
    std::map<ZZ, float> entropies;

    int tot = EntryVector.length();
    for (auto El : EntryVector)
    {
        if(entropies.count(El._ZZ_p__rep)){
            entropies[El._ZZ_p__rep] += 1;
        }
        else{
            entropies[El._ZZ_p__rep] = 1;
        }
    }



    // for (auto it = entropies.begin(); it != entropies.end(); it++)
    // {
    //    // auto p = it->second/tot;
    //   //  entropy -= (p * log2(p));
    // }


    return entropies;
}



void PopulateRandomMatrix(Mat<ZZ_p> &gfMatrix, Vec<Vec<float>>& probabilities){
    int mod = conv<int>(ZZ_p::modulus());


    if(probabilities[0].length() != mod) {

        throw;
    }

    Vec<Vec<float>> distrib {};



    for (size_t i = 0; i < probabilities.length(); i++)
    {
        auto new_vec = Vec<float>{};
        new_vec.append(probabilities[i][0]);
        distrib.append(new_vec);

    }

    for (size_t j = 0; j< probabilities.length(); j++)
    {
        for (size_t i = 1; i < mod; i++)
        {
            distrib[j].append(distrib[j][i-1] + probabilities[j][i]);
            /* code */
        }
    }

    for (size_t i = 0; i < gfMatrix.NumRows(); i++)
    {
        for (size_t j = 0; j < gfMatrix.NumCols(); j++)
        {
            ZZ_p v = GetRandomVWithP(distrib[i]);

            gfMatrix.put(i, j, v);
        }
        auto freq = GetFrequency(gfMatrix[i]);
        for (auto f : freq){
            std::cout << f.first << "_" << f.second << std::endl;

        }
    }



}

float CalculateEntropy(Vec<ZZ_p>& EntryVector){
    float entropy=0;
    std::map<ZZ, float> entropies;

    int tot = EntryVector.length();
    for (auto El : EntryVector)
    {
        if(entropies.count(El._ZZ_p__rep)){
            entropies[El._ZZ_p__rep] += 1;
        }
        else{
            entropies[El._ZZ_p__rep] = 1;
        }
    }



    for (auto it = entropies.begin(); it != entropies.end(); it++)
    {
        auto p = it->second/tot;
        entropy -= (p * log2(p));
    }


    return entropy;


}
float CalculateEntropy(std::vector<int>& EntryVector){
    float entropy=0;
    std::map<int, float> entropies;

    int tot = EntryVector.size();
    for (auto El : EntryVector)
    {
        if(entropies.count(El)){
            entropies[El] += 1;
        }
        else{
            entropies[El] = 1;
        }
    }



    for (auto it = entropies.begin(); it != entropies.end(); it++)
    {
        auto p = it->second/tot;
        entropy -= (p * log2(p));
    }


    return entropy;


}

void CalculateEntropyForAll(std::map<float, Mat<ZZ_p>>& entropies,Mat<ZZ_p>& X , int P, int K, Vec<ZZ_p>* indexes = nullptr){
    auto allocated = false;
    if(indexes == nullptr){
        indexes = new Vec<ZZ_p>{};
        allocated = true;
    }

    if(K == 0){
        if(std::count(indexes->begin(), indexes->end(), ZZ_p{0}) == indexes->length()){
            return;
        }
        Mat<ZZ_p> Row = Mat<ZZ_p>(INIT_SIZE, 1, indexes->length());
        Row[0] = *indexes;
        auto e = CalculateEntropy((Row * X)[0]);
        std::cout << Row << " " << e << std::endl;
        entropies[e] = Row;
        
        return;
    }

    for (size_t i = 0; i < P; i++)
    {
        auto index_to_send = Vec<ZZ_p>{*indexes};
        index_to_send.append(ZZ_p{i});
        CalculateEntropyForAll(entropies, X, P, K-1, &index_to_send);
    }

    if(allocated){
        delete indexes;
    }
    
}


int main()
{
    const int P = 2; //Galois order
    const int K = 6; //Number of sources
    const int T = 10000; // Observation length
    const int pk = pow(P, K);

   // auto probs = rand() % K + P;
    auto p = ZZ_p{1};
    auto z = ZZ_p::zero();

    ZZ_p::init(ZZ(P));
    Mat<ZZ_p> gfMatrix = Mat<ZZ_p>(INIT_SIZE , K, T);

    Vec<float> source_1 {};
    source_1.append(0.2f);
    source_1.append(0.8f);


    Vec<float> source_2 {};
    source_2.append(0.9f);
    source_2.append(0.1f);

    Vec<float> source_3 {};
    source_3.append(0.7f);
    source_3.append(0.3f);

    Vec<float> source_4 {};
    source_4.append(0.4f);
    source_4.append(0.6f);

    Vec<float> source_5 {};
    source_5.append(0.2f);
    source_5.append(0.8f);

    Vec<float> source_6 {};
    source_6.append(0.6f);
    source_6.append(0.4f);


    Vec<Vec<float>> sources;

    sources.append(source_1);
    sources.append(source_2);
    sources.append(source_3);
    sources.append(source_4);
    sources.append(source_5);
    sources.append(source_6);



    PopulateRandomMatrix(gfMatrix, sources);


    Mat<ZZ_p> mixMatrix = Mat<ZZ_p>(INIT_SIZE, K, K);
    PopulateRandomMatrix(mixMatrix);


    auto triu = mixMatrix;
    auto tril = mixMatrix;
    for (size_t i = 0; i < mixMatrix.NumRows(); i++)
    {

        for (size_t j = 0; j < i; j++)
        {
            triu.put(i, j, z);
        }
        if(triu.get(i, i) == 0){
            triu.put(i,i, p);
        }

        for (size_t j = i; j < K; j++)
        {
            tril.put(i, j, z);

        }

        tril.put(i, i, p);


    }

    auto A = triu * tril;
    auto X = A *gfMatrix;

    Mat<ZZ_p> FinalMatrix(INIT_SIZE, K, K);
    std::map<float, Mat<ZZ_p>> entropies;


    CalculateEntropyForAll(entropies, X, P, K);

    for(auto c: entropies){
        std::cout << c.second[0] << " " << c.first << " || "  << std::endl;
    }


    auto it = entropies.cbegin();
    for (size_t i = 0; i < K; i++)
    {
        while(true){
            
            std::cout << it->second << " " << i << std::endl;

            for (size_t z = 0; z < it->second[0].length(); z++)
            {
                FinalMatrix.put(i, z, it->second[0][z]);
            }
            auto m = FinalMatrix;
            
            std::cout << FinalMatrix << std::endl;

            auto r = gauss(m); // Calculo do rank
            std::cout << r << std::endl;

           // std::cout <<  it->first << std::endl;
            it = std::next(it, 1);
            if(r == long(i + 1)){
                break;
            }
            else{
                if(it == entropies.end()){
                    break;
                }
            }
        }
    }

    std::cout << "Final matrix:" << std::endl << FinalMatrix << std::endl;
    std::cout << "A matrix:" << std::endl << A << std::endl;
    std::cout << FinalMatrix * A << std::endl;

   // int x=  conv<int>(p);






















   // random(&Matrix);

}