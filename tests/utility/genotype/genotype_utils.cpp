#include <biovoltron/utility/genotype/genotype_utils.hpp>
#include <iostream> //debug
#include <catch.hpp>

TEST_CASE("Genotype utils") {

    SECTION("get_vcf_genotypes") {
        for(int i=0;i<=biovoltron::GenotypeUtils::MAX_ALLELE_COUNT;i++){
            const auto genotypes = biovoltron::GenotypeUtils::get_vcf_genotypes(i);
            REQUIRE(genotypes.size() == (i * (i + 1)) / 2);
            for(int j=0;j<genotypes.size();j++){
                const auto& genotype = genotypes[j];
                REQUIRE(genotype.first <= genotype.second);
                REQUIRE(genotype.first >= 0);
                REQUIRE(genotype.second < i);
            }
            for(int i=0;i<genotypes.size();i++){
                auto [j,k]=genotypes[i];
                REQUIRE(i==k*(k+1)/2+j);
            }
            REQUIRE(std::is_sorted(genotypes.begin(),genotypes.end(),[](const auto& a, const auto& b) {
                return std::tie(a.second, a.first) < std::tie(b.second, b.first);
            }));
        }
    };

    SECTION("get_raw_genotypes") {
        for(int i=0;i<=biovoltron::GenotypeUtils::MAX_ALLELE_COUNT;i++){
            const auto genotypes = biovoltron::GenotypeUtils::get_raw_genotypes(i);
            REQUIRE(genotypes.size() == (i * (i + 1)) / 2);
            for(int j=0;j<genotypes.size();j++){
                const auto& genotype = genotypes[j];
                REQUIRE(genotype.first <= genotype.second);
                REQUIRE(genotype.first >= 0);
                REQUIRE(genotype.second < i);
            }
            REQUIRE(std::is_sorted(genotypes.begin(),genotypes.end()));
        }
    };

    SECTION("get_genotype_size") {
        for(int i=0;i<=biovoltron::GenotypeUtils::MAX_ALLELE_COUNT;i++){
            const auto size = biovoltron::GenotypeUtils::get_genotype_size(i);
            REQUIRE(size == (i * (i + 1)) / 2);
        }
    };

    SECTION("get_allele_size") {
        for(int i=0;i<=biovoltron::GenotypeUtils::MAX_ALLELE_COUNT;i++){
            const auto size = biovoltron::GenotypeUtils::get_allele_size((i * (i + 1)) / 2);
            REQUIRE(size == i);
        }
    };

    SECTION("to_vcf_order") {
       for(int i=0;i<=biovoltron::GenotypeUtils::MAX_ALLELE_COUNT;i++){
            const auto size = biovoltron::GenotypeUtils::get_allele_size((i * (i + 1)) / 2);
            REQUIRE(size == i);
            std::vector<double> pls((i * (i + 1)) / 2, 0.0);
            for(int ii=0;ii<pls.size();ii++)pls[ii]=ii;
            const auto vcf_pls = biovoltron::GenotypeUtils::to_vcf_order(pls);
            REQUIRE(vcf_pls.size() == pls.size());
            for(int j=0;j<vcf_pls.size();j++){
                int v=(int)std::round(vcf_pls[j]);
                int y;
                for(y=0;y*(y+1)<=2*j;y++);
                --y;
                int x=j-y*(y+1)/2;
               int v2=(size+size-x+1)*x/2+y-x;
               REQUIRE(v==v2);
            }
        }
    };

    SECTION("gls_to_pls") {
        std::vector<double> gls = {0.0, -1.0, -2.0, -3.0, -4.0};
        const auto pls = biovoltron::GenotypeUtils::gls_to_pls(gls);
        REQUIRE(pls.size() == gls.size());
        for (size_t i = 0; i < pls.size(); ++i) {
            REQUIRE(pls[i] == std::round(-10 * (gls[i] - std::ranges::max(gls))));
        }
    }
}