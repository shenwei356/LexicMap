// Original float64 kernels retained only for differential tests and benchmarks.
// Keep their scoring, packing and tie-breaking independent of the optimized kernels.
package cmd

func referenceChaining2(ce *Chainer2, subs *[]*SubstrPair) (*[]*Chain2Result, int, int, int, int, int, int, int) {
	n := len(*subs)
	if n == 1 { // for one seed, just check the seed weight
		sub := (*subs)[0]
		slen := int(sub.Len)
		if slen >= ce.options.MinScore && slen >= ce.options.MinAlignLen { // the length of anchor (max 32)
			paths := poolChains2.Get().(*[]*Chain2Result)
			path := poolChain2.Get().(*Chain2Result)
			path.Reset()
			qe := int(sub.QBegin) + slen - 1           // end
			te := int(sub.TBegin) + slen - 1           // end
			qb, tb := int(sub.QBegin), int(sub.TBegin) // in case there's only one anchor
			path.QBegin, path.QEnd = qb, qe
			path.TBegin, path.TEnd = tb, te
			path.MatchedBases = slen
			path.PIdent = 100
			path.AlignedBasesQ = slen
			path.NAnchors++
			*paths = append(*paths, path)
			return paths, slen, slen, slen, qb, qe, tb, te
		}
		return nil, 0, 0, 0, 0, 0, 0, 0
	}
	var i, j int
	bandBase := int32(ce.options.BandBase) // band size of banded-DP
	bandCount := ce.options.BandCount
	var _bCount int
	var _bBase int32
	maxscoresIdxs := &ce.maxscoresIdxs
	*maxscoresIdxs = (*maxscoresIdxs)[:0]
	*maxscoresIdxs = append(*maxscoresIdxs, uint64((*subs)[0].Len)<<32)
	var s, m, M, g float64
	var mj, Mi int
	var a, b *SubstrPair
	maxGap := float64(ce.options.MaxGap)
	var aQBegin, aTBegin, bQBegin, bTBegin int32
	var qDiff, tDiff int32
	for i = 1; i < n; i++ {
		a = (*subs)[i] // current seed/anchor
		m, mj = float64(a.Len), i
		aQBegin, aTBegin = a.QBegin, a.TBegin
		j = i
		_bCount = 0
		for {
			j--
			if j < 0 {
				break
			}
			b = (*subs)[j] // previous seed/anchor
			bQBegin, bTBegin = b.QBegin, b.TBegin
			if bQBegin == aQBegin || bTBegin > aTBegin {
				continue
			}
			_bCount++
			_bBase = aQBegin - bQBegin - int32(b.Len)
			if !(_bBase <= bandBase || _bCount <= bandCount) {
				break
			}
			qDiff = aQBegin - bQBegin
			tDiff = aTBegin - bTBegin
			if qDiff < 0 {
				qDiff = -qDiff
			}
			if tDiff < 0 {
				tDiff = -tDiff
			}
			if qDiff > tDiff {
				g = float64(qDiff - tDiff)
			} else {
				g = float64(tDiff - qDiff)
			}
			if g > maxGap { // limit the gap. necessary?
				continue
			}
			s = float64((*maxscoresIdxs)[j]>>32) + float64(b.Len) - g // compute the score
			if s >= m {                                               // update the max score of current seed/anchor
				m = s
				mj = j
			}
		}
		*maxscoresIdxs = append(*maxscoresIdxs, uint64(m)<<32|uint64(mj))
		if m > M { // the biggest score in the whole score matrix
			M, Mi = m, i
		}
	}
	minScore := float64(ce.options.MinScore)
	minAlignLen := ce.options.MinAlignLen
	if M < minScore {
		return nil, 0, 0, 0, 0, 0, 0, 0
	}
	paths := poolChains2.Get().(*[]*Chain2Result)
	var nMatchedBases, nAlignedBasesQ, nAlignedBasesT int
	_, qB, qE, tB, tE := chainARegion(
		subs,
		maxscoresIdxs,
		0,
		minScore,
		minAlignLen,
		ce.options.MinIdentity,
		paths,
		&nMatchedBases,
		&nAlignedBasesQ,
		&nAlignedBasesT,
		Mi,
		nil, // &ce.bounds,
		ce.options.HeuristicKmerPidentThreshold,
	)
	if len(*paths) == 0 {
		recycleChaining2ResultSlice(paths)
		return nil, 0, 0, 0, 0, 0, 0, 0
	}
	return paths, nMatchedBases, nAlignedBasesQ, nAlignedBasesT, qB, qE, tB, tE
}

func referenceChaining3(ce *Chainer3, subs *[]*SubstrPair) *Chain3Result {
	n := len(*subs)
	var i, j int
	bandBase := int32(ce.options.BandBase) // band size of banded-DP
	bandCount := ce.options.BandCount
	var _bCount int
	var _bBase int32
	maxscoresIdxs := &ce.maxscoresIdxs
	*maxscoresIdxs = (*maxscoresIdxs)[:0]
	var s, m, M, g, d float64
	var mj, Mi int
	var a, b *SubstrPair
	maxGap := float64(ce.options.MaxGap)
	maxDistance := float64(ce.options.MaxDistance)
	a = (*subs)[0]
	m = float64(a.Len) - distance2(sub0, a) - gap2(sub0, a)
	*maxscoresIdxs = append(*maxscoresIdxs, int64(m)<<32)
	for i = 1; i < n; i++ {
		a = (*subs)[i] // current seed/anchor
		m, mj = float64(a.Len)-distance2(sub0, a)-gap2(sub0, a), i
		j = i
		_bCount = 0
		for {
			j--
			if j < 0 {
				break
			}
			b = (*subs)[j] // previous seed/anchor
			if b.QBegin == a.QBegin || b.TBegin > a.TBegin {
				continue
			}
			_bCount++
			_bBase = a.QBegin - b.QBegin - int32(b.Len)
			if !(_bBase <= bandBase || _bCount <= bandCount) {
				break
			}
			d = distance2(a, b)
			if d > maxDistance {
				continue
			}
			g = gap2(a, b)
			if g > maxGap {
				continue
			}
			s = float64((*maxscoresIdxs)[j]>>32) + float64(b.Len) - d - g // compute the score
			if s >= m {                                                   // update the max score of current seed/anchor
				m = s
				mj = j
			}
		}
		*maxscoresIdxs = append(*maxscoresIdxs, int64(m)<<32|int64(mj))
		if m > M { // the biggest score in the whole score matrix
			M, Mi = m, i
		}
	}
	minScore := float64(ce.options.MinScore)
	minAlignLen := ce.options.MinAlignLen
	if M < minScore {
		return nil
	}
	var nMatchedBases, nAlignedBasesQ, nAlignedBasesT int
	i = Mi
	var qb, qe, tb, te int32 // the bound (0-based)
	var sub *SubstrPair
	var beginOfNextAnchor int
	var pident float64
	firstAnchorOfAChain := true
	var nAnchors int
	for {
		j = int((*maxscoresIdxs)[i] & 4294967295) // previous seed
		if j < 0 {                                // the first anchor is not in current region
			break
		}
		sub = (*subs)[i]
		nAnchors++
		if firstAnchorOfAChain {
			firstAnchorOfAChain = false
			qe = int32(sub.QBegin) + int32(sub.Len) - 1   // end
			te = int32(sub.TBegin) + int32(sub.Len) - 1   // end
			qb, tb = int32(sub.QBegin), int32(sub.TBegin) // in case there's only one anchor
			nMatchedBases += int(sub.Len)
		} else {
			qb, tb = int32(sub.QBegin), int32(sub.TBegin) // begin
			if int(sub.QBegin)+int(sub.Len)-1 >= beginOfNextAnchor {
				nMatchedBases += beginOfNextAnchor - int(sub.QBegin)
			} else {
				nMatchedBases += int(sub.Len)
			}
		}
		beginOfNextAnchor = int(sub.QBegin)
		if i == j { // the path starts here
			if firstAnchorOfAChain { // sadly, there's no anchor added.
				break
			}
			nAlignedBasesQ += int(qe) - int(qb) + 1
			if nAlignedBasesQ < minAlignLen {
				firstAnchorOfAChain = true
				break
			}
			nAlignedBasesT += int(te) - int(tb) + 1
			pident = float64(nMatchedBases) / float64(max(nAlignedBasesQ, nAlignedBasesT)) * 100
			if pident < 15 {
				firstAnchorOfAChain = true
				break
			}
			if pident > 100 {
				pident = 100
			}
			path := poolChain3.Get().(*Chain3Result)
			path.Reset()
			path.AlignedBasesQ = nAlignedBasesQ
			path.AlignedBasesT = nAlignedBasesT
			path.MatchedBases = nMatchedBases
			path.PIdent = pident
			path.QBegin, path.QEnd = int(qb), int(qe)
			path.TBegin, path.TEnd = int(tb), int(te)
			firstAnchorOfAChain = true
			return path
		}
		i = j
	}
	return nil
}
