package cmd

import (
	"slices"
	"testing"
)

// chunkTestDetail builds a minimal SimilarityDetail with the fields the
// merge/finalize/compare helpers read: SimilarityScore, SeqID, and one chain
// covering the given query/subject regions.
func chunkTestDetail(score float64, seqid string, qb, qe, tb int) *SimilarityDetail {
	chains := []*Chain2Result{{QBegin: qb, QEnd: qe, TBegin: tb, TEnd: tb + 100}}
	return &SimilarityDetail{
		SimilarityScore: score,
		SeqID:           []byte(seqid),
		Similarity:      &SeqComparatorResult{Chains: &chains},
	}
}

func chunkTestResult(batch uint64, sds ...*SimilarityDetail) *SearchResult {
	list := sds
	return &SearchResult{BatchGenomeIndex: batch, SimilarityDetails: &list}
}

// A chunk group must report completion only after every dispatched chunk
// signals done — never earlier, and results accumulate in arrival order.
func TestPendingGenomeChunksCompletion(t *testing.T) {
	p := &pendingGenomeChunks{expected: 3}
	p.addResult(chunkTestResult(1, chunkTestDetail(5, "a", 0, 10, 0)))
	if p.markDone() {
		t.Fatal("complete after 1 of 3 done signals")
	}
	p.addResult(chunkTestResult(2, chunkTestDetail(5, "b", 0, 10, 0)))
	if p.markDone() {
		t.Fatal("complete after 2 of 3 done signals")
	}
	// the third chunk may have produced no result, but still signals done
	if !p.markDone() {
		t.Fatal("not complete after 3 of 3 done signals")
	}
	if len(p.results) != 2 {
		t.Fatalf("results = %d, want 2", len(p.results))
	}
}

// mergeChunkResults must concatenate all chunks' SimilarityDetails into the
// first result (in argument order) and leave it usable by finalize.
func TestMergeChunkResults(t *testing.T) {
	idx := &Index{}
	// the later-arriving chunk has the smaller BatchGenomeIndex; the merged
	// result must still keep it as the canonical genome index
	a := chunkTestResult(2, chunkTestDetail(9, "s2", 0, 50, 0))
	b := chunkTestResult(1, chunkTestDetail(8, "s0", 50, 80, 0), chunkTestDetail(7, "s1", 80, 100, 0))

	got := idx.mergeChunkResults([]*SearchResult{a, b})
	if got != a {
		t.Fatal("merged result is not the first chunk's result")
	}
	if got.BatchGenomeIndex != 1 {
		t.Fatalf("merged BatchGenomeIndex = %d, want min 1", got.BatchGenomeIndex)
	}
	if len(*got.SimilarityDetails) != 3 {
		t.Fatalf("merged SimilarityDetails = %d, want 3", len(*got.SimilarityDetails))
	}
	ids := []string{string((*got.SimilarityDetails)[0].SeqID),
		string((*got.SimilarityDetails)[1].SeqID),
		string((*got.SimilarityDetails)[2].SeqID)}
	if !slices.Equal(ids, []string{"s2", "s0", "s1"}) {
		t.Fatalf("merge order = %v, want [s2 s0 s1]", ids)
	}
	if b.SimilarityDetails != nil {
		t.Fatal("donor chunk's details were not recycled")
	}
}

// finalizeGenomeResult must compute AlignedFraction from the union of query
// regions across all (possibly merged) chunk alignments — this is the
// coverage the per-genome filter applies on chunked indexes.
func TestFinalizeGenomeResultCoverage(t *testing.T) {
	idx := &Index{opt: &IndexSearchingOptions{MinQueryAlignedFractionInAGenome: 50}}

	// two chains on different seqids covering [0,39] and [50,89] of a 100bp
	// query (QEnd inclusive) -> 80% genome-level coverage even though each
	// HSP covers less
	r := chunkTestResult(1,
		chunkTestDetail(9, "s0", 0, 39, 10),
		chunkTestDetail(8, "s1", 50, 89, 20))
	if !idx.finalizeGenomeResult(r, 100) {
		t.Fatal("result filtered out although coverage is 80% >= 50%")
	}
	if r.AlignedFraction != 80 {
		t.Fatalf("AlignedFraction = %v, want 80", r.AlignedFraction)
	}

	// a result below the threshold is filtered out (and recycled)
	r2 := chunkTestResult(1, chunkTestDetail(9, "s0", 0, 29, 10))
	if idx.finalizeGenomeResult(r2, 100) {
		t.Fatal("result kept although coverage is 30% < 50%")
	}
}

// finalizeGenomeResult must also sort the (merged) alignments canonically:
// score desc, then seqid, then subject position — independent of the order
// chunk results arrived in.
func TestFinalizeGenomeResultSorts(t *testing.T) {
	idx := &Index{opt: &IndexSearchingOptions{MinQueryAlignedFractionInAGenome: 0}}

	r := chunkTestResult(1,
		chunkTestDetail(9, "s1", 0, 90, 50),
		chunkTestDetail(9, "s0", 0, 90, 10),
		chunkTestDetail(8, "s0", 0, 80, 30))
	if !idx.finalizeGenomeResult(r, 100) {
		t.Fatal("unexpectedly filtered")
	}
	got := make([][2]interface{}, 3)
	for i, sd := range *r.SimilarityDetails {
		got[i] = [2]interface{}{string(sd.SeqID), sdTBegin(sd)}
	}
	want := [][2]interface{}{{"s0", 10}, {"s1", 50}, {"s0", 30}}
	if !slices.Equal(got, want) {
		t.Fatalf("sd order = %v, want %v", got, want)
	}
}

// compareSimilarityDetails: primary key score desc; ties broken by seqid
// then subject position, so both buffered and streamed output order the
// same way regardless of completion order.
func TestCompareSimilarityDetails(t *testing.T) {
	hi := chunkTestDetail(10, "z", 0, 10, 0)
	lo := chunkTestDetail(1, "a", 0, 10, 0)
	if compareSimilarityDetails(hi, lo) >= 0 {
		t.Fatal("higher score did not sort first")
	}
	a := chunkTestDetail(5, "a", 0, 10, 0)
	b := chunkTestDetail(5, "b", 0, 10, 0)
	if compareSimilarityDetails(a, b) >= 0 {
		t.Fatal("score tie not broken by seqid")
	}
	near := chunkTestDetail(5, "a", 0, 10, 10)
	far := chunkTestDetail(5, "a", 0, 10, 90)
	if compareSimilarityDetails(near, far) >= 0 {
		t.Fatal("score+seqid tie not broken by subject position")
	}
}
