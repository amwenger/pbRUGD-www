from __future__ import division
import collections
import functools
import itertools
import math

HpoTerm = collections.namedtuple("HpoTerm", ["id","name","definition"])
HpoAnnotation = collections.namedtuple("HpoAnnotation", ["term","gene","conditions"])

class Phenotyper:
    @staticmethod
    def term_closure(term_list, term_to_relatives):
        """Take the closure over a list of terms.  The closure is defined as the
        terms in the list and all its ancestors or descendants depending on the dag passed."""
        closure = set()

        def helper(term):
            # avoid retraversing a subtree
            if term not in closure:
                closure.add(term)
                for parent in term_to_relatives.get(term, set()):
                    helper(parent)

        # Traverse starting with each term in the term list
        for term in term_list:
            helper(term)
        return closure

    @staticmethod
    def canonicalize_annotations(term_to_parents, gene_to_annotations, term_to_annotations):
        # Take the "closure" over gene<->term annotations such that
        # if a gene is annotated with term T, it is also annotated with
        # all the ancestors of T up to the root of the DAG.
        for gene_id in gene_to_annotations:
            direct_terms = set([a.term for a in gene_to_annotations[gene_id].values()])
            closure = Phenotyper.term_closure(direct_terms, term_to_parents)

            # Add annotations for terms that are in the closure but not in the
            # directly annotated terms.
            for term_id in closure.difference(direct_terms):
                annotation = HpoAnnotation(term_id, gene_id, [])
                gene_to_annotations[gene_id][term_id] = annotation
                term_to_annotations[term_id][gene_id] = annotation
        return gene_to_annotations, term_to_annotations

    @staticmethod
    def compute_term_info_content(hpoterms, term_to_parents, gene_to_annotations, term_to_annotations):
        # Calculate the information content of each term.
        # We define two concepts of information content:
        #    1. "raw": the standard meaning -log2(<genes with term>/<annotated genes>)
        #    2. "marginal": the information content that a term contributes
        #       beyond the information content of its parent terms
        annotatedGeneCt = len(gene_to_annotations)
        termRawIc = dict()
        for term in hpoterms.values():
            annotations = term_to_annotations.get(term.id,[])
            # If term is not annotated, its raw ic is 0 (or Infinity, but that would mess things up so set to 0)
            termRawIc[term.id] = -math.log(len(annotations)/annotatedGeneCt)/math.log(2) if len(annotations) else 0

        # Define marginal information content as the raw information content of a term
        # minus the raw information content of a meta-term that has a gene list that is
        # the intersection of the gene lists of the parents of the term.  The intersection
        # of the gene lists for the parent terms is guaranteed to be at least as large as
        # the gene list of the term itself because of closure.
        termMarginalIc = dict()
        for term in hpoterms.values():
            if term.id not in term_to_annotations:
                termMarginalIc[term.id] = 0
            else:
                parentTerms = term_to_parents.get(term.id, [])
                if len(parentTerms) == 0:
                    parentIc = 0
                else:
                    parentGeneSets = [set([a.gene for g,a in term_to_annotations[t].items()]) for t in parentTerms]
                    parentIntersection = functools.reduce(lambda a,b: a.intersection(b), parentGeneSets)
                    parentIc = -math.log(len(parentIntersection)/annotatedGeneCt)/math.log(2)
                termMarginalIc[term.id] = termRawIc[term.id] - parentIc

        return termMarginalIc

    def __init__ (self, termsfile, dagfile, annotationsfile):
        """Initialize the phenotyper with HPO terms, DAG, and gene<->term annotations."""

        # terms
        self._hpoTerms = {}
        f = open(termsfile)
        for l in f:
            (term_id, term_name, term_defn) = l.rstrip("\n").split("\t")
            self._hpoTerms[term_id] = HpoTerm(term_id,term_name,term_defn)
        f.close()

        # bottom-up DAG (child to parents) and top-down DAGs (parent to children)
        self._buHpoDag = collections.defaultdict(list)
        self._tdHpoDag = collections.defaultdict(list)
        f = open(dagfile)
        for l in f:
            (child_id, parent_id) = l.rstrip("\n").split("\t")
            self._buHpoDag[child_id].append(parent_id)
            self._tdHpoDag[parent_id].append(child_id)

        # gene<->terms annotations and term<->genes annotations
        self._hpoGeneToAnnotations = collections.defaultdict(dict)
        self._hpoTermToAnnotations = collections.defaultdict(dict)
        f = open(annotationsfile)
        for l in f:
            (gene_id, term_id, conditions) = l.rstrip("\n").split("\t")
            annotation = HpoAnnotation(term_id, gene_id, conditions.split(","))
            self._hpoGeneToAnnotations[gene_id][term_id] = annotation
            self._hpoTermToAnnotations[term_id][gene_id] = annotation
        f.close()
        Phenotyper.canonicalize_annotations(self._buHpoDag, self._hpoGeneToAnnotations,self._hpoTermToAnnotations)

        self._termMarginalIc = Phenotyper.compute_term_info_content(self._hpoTerms, self._buHpoDag,
                                                                    self._hpoGeneToAnnotations, self._hpoTermToAnnotations)

        # Define body system indexes
        self._ixToSystem = [s for s,t in bodySystemTerms.items()]
        self._otherIx = len(self._ixToSystem) - 1

        # Memoize _termToSystemIx for efficiency
        self._mzTermToSystemIx = dict()
        # Start by adding base cases.
        for ix,(system,terms) in enumerate(bodySystemTerms.items()):
            for term in terms:
                self._mzTermToSystemIx[term] = ix

    def term_ancestor_closure(self, term_list):
        """Take the closure over a list of terms.  The closure is defined as the
        terms in the list and all ancestors up to the root of the DAG."""
        return Phenotyper.term_closure(term_list, self._buHpoDag)

    def term_descendant_closure(self, term_list):
        """Identify all of the terms in term list or its descendants down
        to the leaves in the DAG."""
        return Phenotyper.term_closure(term_list, self._tdHpoDag)

    def _phenotype_information(self, term_list):
        """Find the information content of a list of terms, which is the sum
        of the marginal information of the terms in the closure of the term list."""
        return sum([self._termMarginalIc.get(t,0) for t in self.term_ancestor_closure(term_list)])

    def phenotype_score(self, term_list1, term_list2):
        """Find the information content shared by two term lists, which is the
        information content of a term list that is the intersection of the closures of
        the two term lists."""
        return self._phenotype_information(self.term_ancestor_closure(term_list1).intersection(self.term_ancestor_closure(term_list2)))

    def all_gene_scores(self, term_list=None, facet_by_condition=False):
        """Return a map from gene ID to the information content shared by the
        term list and the annotations applied to the gene, for all annotated genes.

        If no term_list is passed, then we return the max score for all the genes
        """
        gene_facet_scores = dict()
        for g,annotations in self._hpoGeneToAnnotations.items():
            gene_facet_terms = dict()
            gene_facet_terms[(g,None)] = set(annotations.keys())

            # define facets
            if facet_by_condition:
                for a in annotations.values():
                    for s in a.conditions:
                        gene_facet_terms.setdefault((g,s), set()).add(a.term)

            # compute scores
            for (facet, facet_terms) in gene_facet_terms.items():
                (g, s) = facet
                gene_facet_scores.setdefault(s, {})[g] = self.phenotype_score(term_list if term_list is not None else facet_terms,
                                                                             facet_terms)

        return gene_facet_scores if facet_by_condition else gene_facet_scores[None]

    def _termToSystemIx(self, hpoTerm):
        """A helper for termToSystem that returns the index of the system
        into which to organize a term."""
        if hpoTerm not in self._mzTermToSystemIx:
            parents = self._buHpoDag.get(hpoTerm, [])
            self._mzTermToSystemIx[hpoTerm] = min([self._otherIx] + [self._termToSystemIx(x) for x in parents])
        return self._mzTermToSystemIx[hpoTerm]


    def termToSystem(self, hpoTerm):
        """Determine the body system into which to organize a term."""
        return self._ixToSystem[self._termToSystemIx(hpoTerm)]


    def hpoTermDict(self, hpoId):
        """Return a dictionary that lists details about an HPO term."""
        return {
            "hpoId": hpoId,
            "name": self._hpoTerms[hpoId].name,
            "defn": self._hpoTerms[hpoId].definition,
            "system": self.termToSystem(hpoId),
            "genesWithTerm": len(self._hpoTermToAnnotations.get(hpoId,[]))
        }


    def geneAnnotations(self, geneId):
        return [a for a in self._hpoGeneToAnnotations.get(geneId,{}).values()]


    def groupBySystem(self, terms):
        """Group a set of terms by body system"""
        result = collections.OrderedDict()
        for bodySystem in BodySystems:
            result[bodySystem] = set()
        for term in terms:
            result[self.termToSystem(term)].add(term)
        return result


    def phenotypeComparison(self, refTerms, compareAnnotations):
        result = dict()
        compareTerms = [a.term for a in compareAnnotations]
        termToConditions = {a.term:a.conditions for a in compareAnnotations}
        compareClosure = self.term_ancestor_closure(compareTerms)
        compareMinimal = self.pruneTerms(compareClosure)

        def systemPhenotypeComparison(refTerms):
            """Helper that compares a single body system group of refTerms with compareTerms."""
            # Identify where the "compare terms" intersect with each refTerm t.
            # Group ref terms based on the intersection.
            sharedToRefTerms = collections.defaultdict(set)
            for t in refTerms:
                shared = frozenset(self.pruneTerms(compareClosure.intersection(self.term_ancestor_closure({t}))))
                sharedToRefTerms[shared].add(t)
            return sharedToRefTerms

        # Group the ref terms by body system.
        rtBySystem = collections.defaultdict(set)
        for rt in refTerms:
            rtBySystem[self.termToSystem(rt)].add(rt)

        # For each compare term that is shared, identify the leaf nodes that are support.
        mzSupport = dict()
        def support(term):
            if term not in mzSupport:
                supportTerms = self.term_descendant_closure({term}).intersection(compareClosure)
                # Identify the minimal terms that provide support;
                # list conditions (disease or animal gene disruption) for any terms, even non-leaf ones.
                mzSupport[term] = { "terms": [self.hpoTermDict(s) for s in self.pruneTerms(supportTerms)],
                    "conditions": list(set(itertools.chain(*[termToConditions[s] for s in supportTerms]))) }
            return mzSupport[term]


        # Within each system, group refTerms by where they intersect with the compare terms.
        for system,systemTerms in rtBySystem.items():
            systemSharedToRefTerms = systemPhenotypeComparison(systemTerms)
            result[system] = list()
            for compareShared, rts in systemSharedToRefTerms.items():
                shared = []
                for t in compareShared:
                    o = self.hpoTermDict(t)
                    o["support"] = support(t)
                    shared.append(o)
                terms = [self.hpoTermDict(t) for t in rts]
                result[system].append({"terms": terms, "shared": shared})

        return result


    def pruneTerms(self, termList):
        """Prune a list of terms to only include "leaf" terms that do not have
        a descendant on the list.  This provides a minimal representation such
        that closure(pruneTerms(termList))==closure(termList)."""
        prune = set() # set of terms to prune
        def helper(term):
            # Mark a term and its ancestors for pruning,
            # if it has not already been pruned.
            if term not in prune:
                prune.add(term)
                for p in self._buHpoDag.get(term,[]):
                    helper(p)

        closure = self.term_ancestor_closure(termList)
        for term in closure:
            for p in self._buHpoDag.get(term,[]):
                helper(p)

        return closure.difference(prune)


# Dictionary that maps from body system to the terms that group within that system.
# The systems are listed in priority order (highest to lowest).  If a term falls within
# multiple systems, it is placed in the highest priority system.
bodySystemTerms = collections.OrderedDict()
bodySystemTerms["Known Disease Pattern"] = ["HP:0000005"]
bodySystemTerms["Growth and development"] = ["HP:0001507", "HP:0001197", "HP:0000240"]
bodySystemTerms["Cancer"] = ["HP:0002664"]
bodySystemTerms["Eyes"] = ["HP:0000478"]
bodySystemTerms["Ears, nose, and throat"] = ["HP:0000598", "HP:0011389", "HP:0011452", "HP:0000366", "HP:0001600", "HP:0000600", "HP:0002087"]
bodySystemTerms["Respiratory"] = ["HP:0002086", "HP:0012252"]
bodySystemTerms["Cardiovascular"] = ["HP:0001626"]
bodySystemTerms["Gastrointestinal"] = ["HP:0011024", "HP:0002019", "HP:0002014", "HP:0011968", "HP:0025031"]
bodySystemTerms["Genitourinary"] = ["HP:0000119", "HP:0000079", "HP:0000078"]
bodySystemTerms["Face and skull"] = ["HP:0000271", "HP:0000929", "HP:0000234"]
bodySystemTerms["Musculoskeletal"] = ["HP:0003011", "HP:0000924", "HP:0003549","HP:0040064"]
bodySystemTerms["Neurologic"] = ["HP:0000707", "HP:0001250", "HP:0011446"]
bodySystemTerms["Skin, hair, and nails"] = ["HP:0001574", "HP:0011138"]
bodySystemTerms["Endocrine"] = ["HP:0000818"]
bodySystemTerms["Metabolism"] = ["HP:0001939"]
bodySystemTerms["Hematologic/lymphatic"] = ["HP:0001871"]
bodySystemTerms["Immune"] = ["HP:0002715"]
bodySystemTerms["Other"] = ["HP:0000769"]

# List of body systems
BodySystems = list(bodySystemTerms.keys())
