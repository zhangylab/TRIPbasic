cd /mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/loop_clusters_250102/250

loopfns=(   /mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops_res250/GM12878_HiChIP_H3K27ac_2Dloops
	        /mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops_res250/K562_HiChIP_H3K27ac_2Dloops \
            /mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops_res250/MyLa_HiChIP_H3K27ac_2Dloops \
			/mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops_res250/HeLa_MicroC_2Dloops \
			/mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops_res250/K562_MicroC_2Dloops \
			/mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops_res250/HeLa_TRIP_H3K27ac_2Dloops)

for bw in 0 1 2 3 4 5 6 7 8 9 10;do
	mkdir -p bw_"$bw"
	cd bw_"$bw"

	#resolution 250
    for loopfn in ${loopfns[@]};do
        Rscript /mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/loop_clusters_250102/identify_loop_clusters.R $loopfn $((bw*250))
    done
	cd ..
done



cd /mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/loop_clusters_250102/500
loopfns=(   /mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops/GM12878_HiChIP_H3K27ac_2Dloops
	        /mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops/K562_HiChIP_H3K27ac_2Dloops \
            /mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops/MyLa_HiChIP_H3K27ac_2Dloops \
			/mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops/HeLa_MicroC_2Dloops \
			/mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops/K562_MicroC_2Dloops \
			/mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/apa_against_loops_241227/loops/HeLa_TRIP_H3K27ac_2Dloops)

for bw in 0 1 2 3 4 5 6 7 8 9 10;do
	mkdir -p bw_"$bw"
	cd bw_"$bw"

	#resolution 500
    for loopfn in ${loopfns[@]};do
        Rscript /mirror/Data/Gengyao_Chen/TRIP_article_figures/apa_benchmarking/loop_clusters_250102/identify_loop_clusters.R $loopfn $((bw*500))
    done
	cd ..
done