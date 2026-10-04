package cmd

import (
	"container/heap"
	"reflect"
	"sync"
	"testing"
)

func TestRecycleSubstrPairsDropsOversizedCapacityAfterDeduplication(t *testing.T) {
	for _, long := range []bool{false, true} {
		limit := thresholdNSubs
		listPool := &sync.Pool{}
		if long {
			limit = thresholdNSubsLong
			listPool = poolSubsLong
		}
		subs := make([]*SubstrPair, 2, limit+1)
		subs[0] = &SubstrPair{QBegin: 1, TBegin: 1, Len: 31}
		subs[1] = &SubstrPair{QBegin: 1, TBegin: 1, Len: 31}
		subPool := &sync.Pool{}
		ClearSubstrPairs(subPool, &subs, 31)
		if len(subs) != 1 {
			t.Fatalf("deduplication retained %d anchors", len(subs))
		}
		RecycleSubstrPairs(subPool, listPool, &subs)
		if subs != nil {
			t.Fatalf("long=%v: shortened oversized seed array was retained", long)
		}
	}
}

func TestClearSubstrPairsClearsDiscardedTail(t *testing.T) {
	keep := &SubstrPair{QBegin: 1, TBegin: 1, Len: 31}
	subs := []*SubstrPair{keep, {QBegin: 1, TBegin: 1, Len: 30}}
	backing := subs
	ClearSubstrPairs(&sync.Pool{}, &subs, 31)
	if len(subs) != 1 || subs[0] != keep || backing[1] != nil {
		t.Fatal("deduplication changed the surviving anchor or retained a discarded anchor")
	}
	RecycleSubstrPairs(&sync.Pool{}, &sync.Pool{}, &subs)
	if backing[0] != nil || backing[1] != nil {
		t.Fatal("recycling retains anchors after deduplication shortened the slice")
	}
}

func TestRecycleTransferredResultSlicesPreservesLiveHSP(t *testing.T) {
	c := &Chain2Result{QBegin: 5, QEnd: 15, QSeq: []byte("ACGT")}
	chains := []*Chain2Result{c, nil, c}
	chainBacking := chains
	chains = chains[:1]
	cr := &SeqComparatorResult{Chains: &chains}
	sd := &SimilarityDetail{Similarity: cr}
	details := []*SimilarityDetail{sd, nil, sd}
	detailBacking := details
	details = details[:1]
	r := &SearchResult{SimilarityDetails: &details}
	results := []*SearchResult{r, nil, r}
	resultBacking := results
	results = results[:1]

	// Each slice has transferred its objects to another owner.
	recycleSearchResultSlice(&results)
	recycleSimilarityDetailSlice(&details)
	recycleChaining2ResultSlice(&chains)
	for i := range chainBacking {
		if chainBacking[i] != nil || detailBacking[i] != nil || resultBacking[i] != nil {
			t.Fatal("recycled slice retains references outside its current length")
		}
	}
	if c.QBegin != 5 || c.QEnd != 15 || string(c.QSeq) != "ACGT" || sd.Similarity != cr {
		t.Fatal("slice recycling modified an HSP still owned by another result")
	}
	// Release the objects once, from their final owner.
	cr.Chains = nil
	RecycleSeqComparatorResult(cr)
	sd.Similarity = nil
	recycleChain2(c)
}

func TestRecycleSearchResultsReleasesNestedReferences(t *testing.T) {
	c := &Chain2Result{
		CIGAR: []byte("1M"),
		QSeq:  make([]byte, 1, maxPooledAlignmentBytes+1),
		TSeq:  []byte("A"), Alignment: []byte("|"),
	}
	chains := []*Chain2Result{c, nil}
	chainBacking := chains
	cr := &SeqComparatorResult{Chains: &chains, TSeq: []byte("A")}
	sd := &SimilarityDetail{Similarity: cr}
	details := []*SimilarityDetail{sd}
	detailBacking := details
	subs := make([]*SubstrPair, 1, thresholdNSubs+1)
	subs[0] = &SubstrPair{Len: 31}
	seedChain := make([]int32, 1, chainerInitSize+1)
	seedChains := make([]*[]int32, 1, thresholdNSubs+1)
	seedChains[0] = &seedChain
	r := &SearchResult{Subs: &subs, Chains: &seedChains, SimilarityDetails: &details}
	results := []*SearchResult{r}
	resultBacking := results
	(&Index{}).RecycleSearchResults(&results)
	if len(results) != 0 || resultBacking[0] != nil || detailBacking[0] != nil || chainBacking[0] != nil {
		t.Fatal("recycled results retain nested references")
	}
	if r.SimilarityDetails != nil || sd.Similarity != nil || cr.Chains != nil || cr.TSeq != nil {
		t.Fatal("recycled objects retain nested results or sequence data")
	}
	if r.Subs != nil || r.Chains != nil || subs != nil || seedChain != nil || seedChains != nil {
		t.Fatal("recycled result retains oversized seed arrays")
	}
	if c.CIGAR != nil || c.QSeq != nil || c.TSeq != nil || c.Alignment != nil {
		t.Fatal("oversized alignment buffers were retained")
	}
}

func TestRecycleChain2KeepsSmallOutputBuffers(t *testing.T) {
	c := &Chain2Result{CIGAR: []byte("1M"), QSeq: []byte("A"), TSeq: []byte("A"), Alignment: []byte("|")}
	capacities := [4]int{cap(c.CIGAR), cap(c.QSeq), cap(c.TSeq), cap(c.Alignment)}
	recycleChain2(c)
	for i, b := range [][]byte{c.CIGAR, c.QSeq, c.TSeq, c.Alignment} {
		if len(b) != 0 || cap(b) != capacities[i] {
			t.Fatal("small output buffer was not reset for reuse")
		}
	}
}

func TestRecycleGenomeQueryAndFragmentsDropsBorrowedReferences(t *testing.T) {
	sequence := []byte("ACGTACGT")
	quality := []byte("IIIIIIII")
	q := &GQuery{seqs: []*[]byte{&sequence}, quals: []*[]byte{&quality}}
	seqBacking, qualBacking := q.seqs, q.quals
	fragments, qualities, _, err := seqs2fragmentsWithQual(&q.seqs, &q.quals, 4, 4)
	if err != nil {
		t.Fatal(err)
	}
	fragmentBacking, qualityBacking := *fragments, *qualities
	recycleFragments(fragments)
	recycleFragments(qualities)
	if string(sequence) != "ACGTACGT" || string(quality) != "IIIIIIII" {
		t.Fatal("fragment recycling modified the source query")
	}
	for _, backing := range [][][]byte{fragmentBacking, qualityBacking} {
		for _, fragment := range backing {
			if fragment != nil {
				t.Fatal("fragment array retains a query buffer")
			}
		}
	}
	RecycleGQuery(q)
	if len(q.seqs) != 0 || len(q.quals) != 0 || seqBacking[0] != nil || qualBacking[0] != nil {
		t.Fatal("query array retains returned sequence or quality buffers")
	}
}

func TestRecycleMergedResultsReleasesRowsAndHeapTail(t *testing.T) {
	record := &SearchResultOfASequence{Sseqid: "contig", Extra: "1M\tA\tA"}
	result := &SearchResultOfAGenome{Query: "q", Qlen: "1", Sgenome: "g", Records: []*SearchResultOfASequence{record}}
	entries := []*SearchResultOfAGenome{result}
	backing := entries
	h := &SearchResultsHeap{entries: &entries}
	popped := heap.Pop(h).(*SearchResultOfAGenome)
	if popped != result || backing[0] != nil || record.Extra == "" {
		t.Fatal("heap pop lost the result or retained its pointer")
	}
	rows := result.Records
	RecycleSearchResultOfAGenome(popped)
	if rows[0] != nil || result.Query != "" || result.Qlen != "" || *record != (SearchResultOfASequence{}) {
		t.Fatal("recycled merge result retains an input row")
	}
	sam := &SearchResultOfASequence2{SEQ: "ACGT", RNAME: "contig", CIGAR: "4M"}
	recycleSearchResultOfASequence2(sam)
	if *sam != (SearchResultOfASequence2{}) {
		t.Fatal("recycled SAM row retains parsed text")
	}
}

func TestSplitClearsUnusedFields(t *testing.T) {
	for _, split := range []func(string, *[]string){
		func(s string, fields *[]string) { stringSplitN(s, "\t", 3, fields) },
		func(s string, fields *[]string) { stringSplitNByByte(s, '\t', 3, fields) },
	} {
		fields := make([]string, 3)
		split("a\tb\tc", &fields)
		backing := fields
		split("d", &fields)
		if !reflect.DeepEqual(fields, []string{"d"}) || backing[1] != "" || backing[2] != "" {
			t.Fatal("short row retains fields from the previous row")
		}
		fields = fields[:3]
		split("e\tf", &fields)
		if !reflect.DeepEqual(fields, []string{"e", "f"}) || backing[2] != "" {
			t.Fatal("reused fields changed the split result")
		}
	}
	fields := make([][]byte, 3)
	bytesSplitN([]byte("a\tb\tc"), []byte("\t"), 3, &fields)
	backing := fields
	bytesSplitN([]byte("d"), []byte("\t"), 3, &fields)
	if len(fields) != 1 || string(fields[0]) != "d" || backing[1] != nil || backing[2] != nil {
		t.Fatal("short byte row retains previous fields")
	}
}
