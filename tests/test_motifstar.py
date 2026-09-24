"""Focused regression tests for MotifSTaR 1.4.7."""
from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("motifstar_1_4_7", ROOT / "motifstar.py")
tool = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = tool
SPEC.loader.exec_module(tool)


def evidence(panel, motifs, row):
    return tool.ReadEvidence(panel, len(motifs), float(row), f"{panel}-{row}",
                             "".join(motifs), "".join(motifs), tuple(motifs))


def infer(rows, sizes, kind="homozygous-size", **options):
    return tool._infer_alleles("sample", kind, list(sizes.values()), rows, sizes,
                              "CAG", set(), set(), tool.RunConfig(**options))


def compound_fixture(vectors, paths=None, counts=(4, 4), locus="FXN", orientation="reference", indel_panel=None):
    motifs = ("A", "GAA", "CCG")[:len(vectors[0])]
    panels = []
    for panel_id, (sizes, count) in enumerate(zip(vectors, counts), 1):
        blocks = [tool.RepeatBlock(index, 0., float(size * len(motif)), float(len(motif)),
                  size, "orange" if index == 0 else "green", 1., len(motif), [], motif * size)
                  for index, (motif, size) in enumerate(zip(motifs, sizes))]
        reads = []
        for row in range(count):
            seqs = [motif * size for motif, size in zip(motifs, sizes)]
            if paths is not None:
                seqs[1] = "".join(paths[panel_id - 1])
            reasons = [""] * len(motifs)
            if panel_id == indel_panel:
                reasons[1] = "target_deletion"
            reads.append(tool.SvgRead(panel_id, float(panel_id*100+row), f"p{panel_id}r{row}",
                                     seqs, reasons, True))
        panels.append(tool.SvgPanel(panel_id, float(panel_id*100), blocks, reads, {}))
    components = [tool.CatalogComponent(i, f"{locus}_{i}", f"chr1:{i*100}-{i*100+30}", motif,
                                       is_primary=i == 1, quantifier="*") for i, motif in enumerate(motifs)]
    catalog = tool.CatalogLocus(locus, locus, "".join(f"({motif})*" for motif in motifs),
                              "GAA", components[1].reference_region, components)
    annotation = tool.Annotation(locus, "+" if orientation == "reference" else "-", orientation, gene=locus)
    return tool.ParsedSvg(panels, []), catalog, annotation


def run_fixture(fixture, mode="combined", selected=()):
    parsed, catalog, annotation = fixture
    root = ROOT / "tmp" / "test_1_4_6_inputs"
    path = root / "sample" / f"sample.{catalog.locus_id}.eh_realigned.reviewer.{catalog.locus_id}.svg"
    context = tool.WorkerContext(str(root), [catalog], [annotation], {}, {}, list(selected),
                                 tool.RunConfig(component_mode=mode))
    tool._worker_initialise(context)
    with mock.patch.object(tool, "parse_reviewer_svg", return_value=parsed):
        result = tool._process_svg_worker(str(path))
    for record in result:
        if record["qc"]["qc_status"] == "ERROR":
            raise AssertionError(record["qc"]["notes"])
    return result


class CompoundVectorTests(unittest.TestCase):
    def test_unreported_companion_prevents_pooling_in_primary_only_fxn(self):
        pure, alt = ("GAA",)*9, ("GAA",)*8+("GAG",)
        record, = run_fixture(compound_fixture([(4, 9), (27, 9)], [pure, alt]))
        self.assertEqual(record["main"]["allele1"], "(GAA)9")
        self.assertEqual(record["main"]["allele2"], "(GAA)8(GAG)1")
        for n in (1, 2):
            self.assertEqual(record["main"][f"consensus_eligible_read_depth_allele{n}"], 4)
            self.assertIn(f"panel_{n}_", record["main"][f"evidence_scope_allele{n}"])

    def test_full_equal_vectors_still_pool(self):
        record, = run_fixture(compound_fixture([(27, 9), (27, 9)]))
        for n in (1, 2):
            self.assertEqual(record["main"][f"allele{n}"], "(GAA)9")
            self.assertEqual(record["main"][f"consensus_eligible_read_depth_allele{n}"], 8)
            self.assertTrue(record["main"][f"evidence_scope_allele{n}"].startswith("pooled_"))

    def test_third_unreported_component_is_checked(self):
        record, = run_fixture(compound_fixture([(27, 9, 2), (27, 9, 3)]))
        self.assertEqual(record["main"]["consensus_eligible_read_depth_allele1"], 4)
        self.assertIn("panel_2_", record["main"]["evidence_scope_allele2"])

    def test_equal_total_sizes_do_not_imply_equal_vectors(self):
        parsed, _, _ = compound_fixture([(4, 27), (27, 4)])
        self.assertEqual(tool._compound_genotype_class(parsed), "heterozygous-size")

    def test_mode_primary_all_and_explicit_component_obey_full_context(self):
        fixture = compound_fixture([(4, 9), (27, 9)], locus="CUSTOM")
        for mode, selected in [("primary", ()), ("all", ()), ("combined", ("CUSTOM_1",))]:
            with self.subTest(mode=mode):
                records = run_fixture(fixture, mode, selected)
                green = next(record for record in records if record["combined"]["locus"] == "CUSTOM_1")
                self.assertEqual(green["main"]["consensus_eligible_read_depth_allele1"], 4)
                self.assertIn("panel_2_", green["main"]["evidence_scope_allele2"])

    def test_combined_reporting_still_appends_catalog_companions_in_panel_order(self):
        record, = run_fixture(compound_fixture([(27, 9), (4, 9)], locus="CUSTOM"))
        self.assertEqual(record["main"]["allele1"], "(A)27(GAA)9")
        self.assertEqual(record["main"]["allele2"], "(A)4(GAA)9")

    def test_two_top_reads_cannot_be_rescued_from_bottom_panel(self):
        record, = run_fixture(compound_fixture([(4, 9), (27, 9)], counts=(2, 5)))
        self.assertEqual(record["main"]["allele1"], "NA")
        self.assertEqual(record["main"]["allele2"], "(GAA)9")
        self.assertNotIn("no_allele_structure_reconstructed", record["qc"]["qc_flags"])

    def test_bottom_indel_is_not_moved_or_rescued(self):
        record, = run_fixture(compound_fixture([(4, 9), (27, 9)], indel_panel=2))
        self.assertEqual(record["main"]["allele1"], "(GAA)9")
        self.assertEqual(record["main"]["allele2"], tool.INDEL_SENTINEL)
        self.assertIn(tool.INDEL_REVIEW_FLAG, record["qc"]["qc_flags"])

    def test_reporting_strand_does_not_swap_panels(self):
        record, = run_fixture(compound_fixture([(4, 17), (27, 9)], orientation="reverse_complement"))
        self.assertEqual(record["main"]["allele1"], "(TTC)17")
        self.assertEqual(record["main"]["allele2"], "(TTC)9")


class PanelOrderTests(unittest.TestCase):
    def test_larger_top_allele_stays_allele_one(self):
        rows = [evidence(1, ("CAG",)*17, i) for i in range(4)]
        rows += [evidence(2, ("CAG",)*9, i) for i in range(3)]
        main, _, flags, _ = infer(rows, {1:17, 2:9}, "heterozygous-size")
        self.assertEqual((main["allele1"], main["allele2"]), ("(CAG)17", "(CAG)9"))
        self.assertEqual((main["read_depth_allele1"], main["read_depth_allele2"]), (4,3))
        self.assertFalse(any("differs_from_eh" in flag for flag in flags))

    def test_same_size_pure_top_interrupted_bottom_not_sorted_by_burden(self):
        rows = [evidence(1, ("CAG",)*3, i) for i in range(4)]
        rows += [evidence(2, ("CAG","CAA","CAG"), i) for i in range(3)]
        main, _, flags, _ = infer(rows, {1:3, 2:3})
        self.assertEqual(main["allele1"], "(CAG)3")
        self.assertEqual(main["allele2"], "(CAG)1(CAA)1(CAG)1")
        self.assertIn("same_size_panel_order_from_read_support", flags)
        self.assertEqual(tool._qc_status(flags), "PASS")

    def test_ambiguous_mixed_panels_keep_distinct_calls_without_review(self):
        rows = [evidence(panel, path, j*10+i) for panel in (1,2)
                for j,path in enumerate([("CAG",)*3, ("CAG","CAA","CAG")]) for i in range(3)]
        main, _, flags, _ = infer(rows, {1:3, 2:3})
        self.assertEqual((main["allele1"],main["allele2"]), ("(CAG)1(CAA)1(CAG)1","(CAG)3"))
        self.assertNotIn("same_size_panel_assignment_unresolved", flags)
        self.assertNotIn("same_size_panel_order_from_read_support", flags)
        self.assertEqual(tool._qc_status(flags), "PASS")
        self.assertNotIn("_panel_", main["evidence_scope_allele1"])
        self.assertEqual((main["read_depth_allele1"],main["read_depth_allele2"]), (6,6))

    def test_partial_pooled_candidate_is_retained_without_order_review(self):
        rows = [evidence(panel, ("CAG",)*3, i) for panel in (1,2) for i in range(3)]
        rows += [evidence(panel, ("CAA","CAG","CAG"), 10) for panel in (1,2)]
        rows += [evidence(panel, ("CAG","CAA","CAG"), 11) for panel in (1,2)]
        rows += [evidence(panel, ("CAG","CAG","CAA"), 12) for panel in (1,2)]
        main, _, flags, notes = infer(rows, {1:3,2:3}, same_size_min_cluster_reads=3,
                                      same_size_max_assignment_distance=0)
        self.assertEqual((main["allele1"],main["allele2"]), ("(CAG)3","NA"))
        self.assertNotIn("same_size_panel_assignment_unresolved", flags)
        self.assertIn("same_size_top_two_coverage_below_threshold", flags)
        self.assertTrue(any("original reporting convention" in note for note in notes))

    def test_identical_pooled_calls_are_unchanged(self):
        rows = [evidence(panel, ("CAG",)*3, i) for panel in (1,2) for i in range(3)]
        main, _, flags, _ = infer(rows, {1:3, 2:3})
        self.assertEqual((main["allele1"],main["allele2"]), ("(CAG)3","(CAG)3"))
        self.assertEqual((main["read_depth_allele1"],main["read_depth_allele2"]), (6,6))

    def test_single_incomplete_candidate_can_belong_to_bottom_panel(self):
        rows = [evidence(1, ("CAA","CAG","CAG"), 1), evidence(1, ("CAG","CAA","CAG"), 2)]
        rows += [evidence(2, ("CAG",)*3, i) for i in range(4)]
        main, _, flags, _ = infer(rows, {1:3, 2:3}, same_size_min_cluster_reads=3,
                                  same_size_max_assignment_distance=0)
        self.assertEqual((main["allele1"],main["allele2"]), ("NA","(CAG)3"))
        self.assertNotIn("no_allele_structure_reconstructed", flags)

    def test_strict_consensus_threshold_is_unchanged(self):
        rows = [evidence(1, path, j*10+i) for j,path in enumerate([("CAG",)*3,("CAG","CAA","CAG")])
                for i in range(3)] + [evidence(2, ("CAG",)*5, i) for i in range(3)]
        main, _, _, _ = infer(rows, {1:3,2:5}, "heterozygous-size")
        self.assertEqual(main["allele1"], "(CAG)3")


if __name__ == "__main__":
    unittest.main(verbosity=2)
