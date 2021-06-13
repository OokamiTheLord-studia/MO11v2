// MO11v2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cmath>
#include <functional>
#include <cassert>
#include <fstream>
#include "net.h"
#include "FillableNet.h"
#include "KMB.h"
#include "CN_Thomas.h"
#include "CN_GaussSiedel.h"
#include "MO11v2.h"

constexpr double t_max{ 1 };
constexpr double b{ 0.1 };
constexpr double d{ 1 };

double analyticSolution(double t, double x)
{
    return 0.5 * exp(((d * t) / (b * b)) - (x / b)) * erfc((((2 * d * t) /b) - x) / (2 * sqrt(d * t)));
}

double maxDiffAtMaxT(MO::FillableNet net1, MO::Net net2)
{
    auto tmax1{ net1.getMatrix()->back() };
    auto tmax2{ net2.getMatrix()->back() };

    assert(tmax1.size() == tmax2.size());

    double maxDiff{ -1 };
    {
        auto itmax1{ tmax1.cbegin() };
        auto itmax2{ tmax2.cbegin() };
        while (itmax1 != tmax1.cend())
        {
            double temp{ abs(*itmax1 - *itmax2) };
            if (temp > maxDiff) maxDiff = temp;

            itmax1++;
            itmax2++;
        }
    }
    return maxDiff;
}

void Zad1(double x_begin, double x_end, const double& t_begin, std::function<double(double)>& start_condition, std::function<double(double)>& edge_condition_derivative_parameter, std::function<double(double)>& edge_condition_function_parameter, std::function<double(double)>& edge_condition_free_function_parameter)
{
    {
        constexpr double h_max{ 0.01 };
        //constexpr double h_max{ 0.095 };
        constexpr double h_min{ 0.1 };
        constexpr double dh{ 0.005 };
        constexpr double GS_h_max{ 0.04 };
        const unsigned int h_count{ static_cast<unsigned int>(std::floor((std::abs(h_min - h_max))/dh)) + 1 };

        std::vector<double> kmb_tmax_err;
        std::vector<double> thom_tmax_err;
        std::vector<double> gs_tmax_err;
        std::vector<double> h_list;
        kmb_tmax_err.reserve(h_count);
        thom_tmax_err.reserve(h_count);
        gs_tmax_err.reserve(h_count);

        for (double h{ h_min }; h >= h_max; h -= dh)
        {
            const double dt1{ h * h };
            const double dt04{ dt1 * 0.4 };

            MO::FillableNet trueSolution(x_begin, x_end, h, t_begin, t_max, dt1, analyticSolution);
            MO::Net KMBSolution(x_begin, x_end, h, t_begin, t_max, dt04, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
            MO::Net ThomasSolution(x_begin, x_end, h, t_begin, t_max, dt1, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
            
            if (h >= GS_h_max)
            {
                MO::Net GaussSiedelSolution(x_begin, x_end, h, t_begin, t_max, dt1, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
                MO::Crank_Nicolson::CN_GaussSiedel GSSolver(1);
                GaussSiedelSolution.solve(&GSSolver);
                gs_tmax_err.push_back(maxDiffAtMaxT(trueSolution, GaussSiedelSolution));
            }

            MO::KMB KMBSolver(0.4);
            MO::Crank_Nicolson::CN_Thomas ThomasSolver(1);
            

            KMBSolution.solve(&KMBSolver);
            ThomasSolution.solve(&ThomasSolver);

            

            h_list.push_back(h);
            kmb_tmax_err.push_back(maxDiffAtMaxT(trueSolution, KMBSolution));
            thom_tmax_err.push_back(maxDiffAtMaxT(trueSolution, ThomasSolution));
            

            std::cout << "Zadanie 1: " << "aktualne h: " << h << std::endl;
        }

        std::ofstream kmbFile;
        std::ofstream thomFile;
        std::ofstream gsFile;
        kmbFile.open("zad1_kmb.csv");
        thomFile.open("zad1_thom.csv");
        gsFile.open("zad1_gs.csv");

        auto ikmb_tmax_err{ kmb_tmax_err.cbegin() };
        auto ithom_tmax_err{ thom_tmax_err.cbegin() };
        auto igs_tmax_err{ gs_tmax_err.cbegin() };
        auto ih_list{ h_list.cbegin() };
        while (ih_list != h_list.cend())
        {
            kmbFile << *ih_list << "," << *ikmb_tmax_err << std::endl;
            thomFile << *ih_list << "," << *ithom_tmax_err << std::endl;
            if (igs_tmax_err != gs_tmax_err.cend())
            {
                gsFile << *ih_list << "," << *igs_tmax_err << std::endl;
                igs_tmax_err++;
            }

            ikmb_tmax_err++;
            ithom_tmax_err++;
            ih_list++;
        }

        kmbFile.close();
        thomFile.close();
        gsFile.close();

    }
}

int main()
{
    double x_end{6 * sqrt(d*t_max)};
    double x_begin{-x_end};
    constexpr double t_begin{ 0 };
    //constexpr double h{ 0.01 };
    //constexpr double dt{ 0.001 };

    std::function<double(double)> temporary_function{ [](double temp) {return 0; } };

    std::function<double(double)> start_condition{ [&](double x) {return x < 0 ? 0 : exp(-x / b); } };
    std::function<double(double)> edge_condition_derivative_parameter{ [](double t) {return 0; } };
    std::function<double(double)> edge_condition_function_parameter{ [](double t) {return 1; } };
    std::function<double(double)> edge_condition_free_function_parameter{ [](double t) {return 0; } };





    /*constexpr double KMBh{ 0.05 };
    constexpr double KMBlambda{ dt / (KMBh * KMBh) };*/

    /*MO::Net KMBSolvedNet(x_begin, x_end, KMBh, t_begin, t_max, dt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
    MO::KMB kmbSolver(KMBlambda);*/
    /*KMBSolvedNet.solve(&kmbSolver);
    KMBSolvedNet.dump("KMB_corrected.csv");*/

    /*constexpr double CNh{ 0.05 };
    constexpr double CNdt{ 0.0025 };*/

    /*MO::Net CNThomasSolvedNet(x_begin, x_end, CNh, t_begin, t_max, CNdt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
    MO::Crank_Nicolson::CN_Thomas CNThomasSolver{ CNdt / (CNh * CNh) };*/
    /*CNThomasSolvedNet.solve(&CNThomasSolver);
    CNThomasSolvedNet.dump("CNThomas.csv");*/

   /* MO::Net CNGaussSiedelSolvedNet(x_begin, x_end, CNh, t_begin, t_max, CNdt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
    MO::Crank_Nicolson::CN_GaussSiedel CNGaussSiedelSolver{ CNdt / (CNh * CNh) };*/
   /* CNGaussSiedelSolvedNet.solve(&CNGaussSiedelSolver);
    CNGaussSiedelSolvedNet.dump("CNGaussSiedel.csv");*/

    //constexpr double testh{ 0.01 };
    //constexpr double testt{ 0.00004 };
    //MO::Net testCNGaussSiedelSolvedNet(x_begin, x_end, testh, t_begin, t_max, testt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
    ///*MO::Crank_Nicolson::CN_Thomas testCNGaussSiedelSolver{ testt / (testh * testh) };*/
    //MO::KMB testKMBSolver(testt / (testh * testh));
    //testCNGaussSiedelSolvedNet.solve(&testKMBSolver);

    //Zad1
    //Zad1(x_begin, x_end, t_begin, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);

    {
        constexpr double KMBh{ 0.01 };
        constexpr double Thomh{ 0.01 };
        constexpr double GSh{ 0.04 };

        constexpr double KMBt{ 0.00004 };
        constexpr double Thomt{ 0.0001 };
        constexpr double GSt{ 0.0016 };

        MO::FillableNet analiticSolutionForKMB(x_begin, x_end, KMBh, t_begin, t_max, KMBt, analyticSolution);
        MO::FillableNet analiticSolutionForThom(x_begin, x_end, Thomh, t_begin, t_max, Thomt, analyticSolution);
        MO::FillableNet analiticSolutionForGS(x_begin, x_end, GSh, t_begin, t_max, GSt, analyticSolution);

        analiticSolutionForKMB.dump("Zad2_AnalKMB.csv");
        analiticSolutionForThom.dump("Zad2_AnalThom.csv");
        analiticSolutionForGS.dump("Zad2_AnalGS.csv");

        return 0;
    }

    //Zad2 i Zad3
    {
        constexpr double KMBh{ 0.01 };
        constexpr double Thomh{ 0.01 };
        constexpr double GSh{ 0.04 };

        constexpr double KMBt{ 0.00004 };
        constexpr double Thomt{ 0.0001 };
        constexpr double GSt{ 0.0016 };

        constexpr double KMBlambda{ KMBt / (KMBh * KMBh) };
        constexpr double Thomlambda{ Thomt / (Thomh * Thomh) };
        constexpr double GSlambda{ GSt / (GSh * GSh) };

        MO::FillableNet analiticSolutionForKMB(x_begin, x_end, KMBh, t_begin, t_max, KMBt, analyticSolution);
        MO::Net KMBSolvedNet(x_begin, x_end, KMBh, t_begin, t_max, KMBt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
        MO::KMB KMBSolver(KMBlambda);

        MO::FillableNet analiticSolutionForThom(x_begin, x_end, Thomh, t_begin, t_max, Thomt, analyticSolution);
        MO::Net ThomSolvedNet(x_begin, x_end, Thomh, t_begin, t_max, Thomt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
        MO::Crank_Nicolson::CN_Thomas ThomasSolver(Thomlambda);

        MO::FillableNet analiticSolutionForGS(x_begin, x_end, GSh, t_begin, t_max, GSt, analyticSolution);
        MO::Net GSSolvedNet(x_begin, x_end, GSh, t_begin, t_max, GSt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
        MO::Crank_Nicolson::CN_GaussSiedel GSSolver(GSlambda);

        KMBSolvedNet.solve(&KMBSolver);
        ThomSolvedNet.solve(&ThomasSolver);
        GSSolvedNet.solve(&GSSolver);

        KMBSolvedNet.dump("Zad2_KMB.csv");
        ThomSolvedNet.dump("Zad2_Thomas.csv");
        GSSolvedNet.dump("Zad2_GS.csv");

        //Zad3
        std::vector<double> KMBerr;
        std::vector<double> Thomerr;
        std::vector<double> GSerr;
        KMBerr.reserve(KMBSolvedNet.getTPositions()->size());
        Thomerr.reserve(ThomSolvedNet.getTPositions()->size());
        GSerr.reserve(GSSolvedNet.getTPositions()->size());

        {
            auto matrixAnal{ analiticSolutionForKMB.getMatrix() };
            auto matrixKMB{ KMBSolvedNet.getMatrix() };
            auto iAnal{ matrixAnal->cbegin() };
            auto iKMB{ matrixKMB->cbegin() };
            while (iKMB != matrixKMB->cend())
            {
                double temp{ -1 };
                auto iElemAnal{ iAnal->cbegin() };
                auto iElemKMB{ iKMB->cbegin() };
                while (iElemKMB != iKMB->cend())
                {
                    double temp2{ std::abs(*iElemAnal - *iElemKMB) };
                    if (temp2 > temp) temp = temp2;

                    iElemAnal++;
                    iElemKMB++;
                }

                KMBerr.push_back(temp);

                iAnal++;
                iKMB++;
            }
        }
        {
            auto matrixAnal{ analiticSolutionForThom.getMatrix() };
            auto matrixThom{ ThomSolvedNet.getMatrix() };
            auto iAnal{ matrixAnal->cbegin() };
            auto iThom{ matrixThom->cbegin() };
            while (iThom != matrixThom->cend())
            {
                double temp{ -1 };
                auto iElemAnal{ iAnal->cbegin() };
                auto iElemThom{ iThom->cbegin() };
                while (iElemThom != iThom->cend())
                {
                    double temp2{ std::abs(*iElemAnal - *iElemThom) };
                    if (temp2 > temp) temp = temp2;

                    iElemAnal++;
                    iElemThom++;
                }

                Thomerr.push_back(temp);

                iAnal++;
                iThom++;
            }
        }
        {
            auto matrixAnal{ analiticSolutionForGS.getMatrix() };
            auto matrixGS{ GSSolvedNet.getMatrix() };
            auto iAnal{ matrixAnal->cbegin() };
            auto iGS{ matrixGS->cbegin() };
            while (iGS != matrixGS->cend())
            {
                double temp{ -1 };
                auto iElemAnal{ iAnal->cbegin() };
                auto iElemGS{ iGS->cbegin() };
                while (iElemGS != iGS->cend())
                {
                    double temp2{ std::abs(*iElemAnal - *iElemGS) };
                    if (temp2 > temp) temp = temp2;

                    iElemAnal++;
                    iElemGS++;
                }

                GSerr.push_back(temp);

                iAnal++;
                iGS++;
            }
        }

        {
            std::ofstream KMBfile;
            KMBfile.open("Zad3_KMB.csv");

            auto times{ KMBSolvedNet.getTPositions() };
            auto itimes{ times->cbegin() };
            auto iKMBerr{ KMBerr.cbegin() };

            while (iKMBerr != KMBerr.cend())
            {
                KMBfile << *itimes << "," << *iKMBerr << std::endl;

                iKMBerr++;
                itimes++;
            }

        }
        {
            std::ofstream Thomfile;
            Thomfile.open("Zad3_Thom.csv");

            auto times{ ThomSolvedNet.getTPositions() };
            auto itimes{ times->cbegin() };
            auto iThomerr{ Thomerr.cbegin() };

            while (iThomerr != Thomerr.cend())
            {
                Thomfile << *itimes << "," << *iThomerr << std::endl;

                iThomerr++;
                itimes++;
            }

        }
        {
            std::ofstream GSfile;
            GSfile.open("Zad3_GS.csv");

            auto times{ GSSolvedNet.getTPositions() };
            auto itimes{ times->cbegin() };
            auto iGSerr{ GSerr.cbegin() };

            while (iGSerr != GSerr.cend())
            {
                GSfile << *itimes << "," << *iGSerr << std::endl;

                iGSerr++;
                itimes++;
            }

        }
        
    }
    
    
}
