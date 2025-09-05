require(InteractionSet)

args = commandArgs(trailingOnly = T)
loopfn = args[1]
bw = as.numeric(args[2])

#1. read loops into GInteractions object
loops = read.table(loopfn, header = F, sep="\t")
an1 = GRanges(seqnames = loops[,1], IRanges(start=loops[,2], end=loops[,3]))
an2 = GRanges(seqnames = loops[,4], IRanges(start=loops[,5], end=loops[,6]))
an1 = resize(an1, width=width(an1)+bw)
an2 = resize(an2, width=width(an2)+bw)
gi = GInteractions(an1, an2)
mcols(gi)$id=1:length(gi)

#2. compute loop overlapping with itself, overlapped loops are considered as clusters
giov = findOverlaps(gi, gi)

#3. overlapped loops in pairs are presented as edge in a undirected graph
require(igraph)
edges = as.data.frame(giov)
nodes = mcols(gi)$id
g = graph_from_data_frame(d = edges, directed = F, vertices = data.frame(name=nodes))
sg = simplify(g)

#4. find the connected components in the graph as clustered loops
comps = components(sg)

#5. connected components are the loop clusters, and component ids are group ids used in bounding box identification which is the smallest rectangle that contains all loops in a cluster
group.by = comps$membership
bb = boundingBox(gi, group.by)

#6. reduce bounding boxes by findOverlaps by considering "within" and "equal" as overlaps
bbov = findOverlaps(bb, bb, type=c("within"))
bbov_reduced1 = bbov[queryHits(bbov) != subjectHits(bbov)]
bbov_to_be_removed1 = queryHits(bbov_reduced1)
bbov = findOverlaps(bb, bb, type=c("equal"))
bbov_reduced2 = bbov[queryHits(bbov) != subjectHits(bbov)]
bbov_to_be_removed2 = queryHits(bbov_reduced2)
bbov_to_be_removed = unique(c(bbov_to_be_removed1, bbov_to_be_removed2))

if(length(bbov_to_be_removed)>0){
    bb = bb[-bbov_to_be_removed]
}
bban1 = anchors(bb, "first")
bban2 = anchors(bb, "second")
bbdf = data.frame( chr1=seqnames(bban1), s1=start(bban1), e1=end(bban1),
                   chr2=seqnames(bban2), s2=start(bban2), e1=end(bban2) )
outfn = paste0(basename(loopfn), ".cluster.bedpe")
write.table(bbdf, outfn, row.names = F, col.names = F, sep="\t", quote = F)