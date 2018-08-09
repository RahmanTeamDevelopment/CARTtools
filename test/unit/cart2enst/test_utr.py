from unittest import TestCase
from mock import MagicMock, patch


class TestUTRSelection(TestCase):


    def _test_utr_selection(self):

        self.assertTrue(1 == 1)


    def test_utr_selection_utr5_single_utr5_exon_number_match(self):

        from selectensts_ import utr

        cart = MagicMock(utr5_exons=[None]*4)
        enst1 = MagicMock(utr5_exons=[None]*2, id_=1)
        enst2 = MagicMock(utr5_exons=[None]*1, id_=2)
        enst3 = MagicMock(utr5_exons=[None]*4, id_=3)
        enst4 = MagicMock(utr5_exons=[None]*3, id_=4)
        enst5 = MagicMock(utr5_exons=[None]*5, id_=5)

        selected, x, y, z = utr._utr_selection_utr5(cart,[enst1, enst2, enst3, enst4, enst5], None)
        self.assertEquals(len(selected), 1)
        self.assertEquals(selected[0].id_, 3)
        self.assertEquals(x, 'single_utr5_exon_number_match')
        self.assertEquals(y, '.')
        self.assertEquals(z, '.')


    @patch('cart2enst_.utr._enst_utr5_encompass_cart_utr5')
    def test_utr_selection_utr5_single_encompassing_utr5(self, mocked_enst_utr5_encompass_cart_utr5):

        from selectensts_ import utr
        mocked_enst_utr5_encompass_cart_utr5.side_effect = lambda cart, enst: enst.id_ == 5

        cart = MagicMock(utr5_exons=[None] * 4)
        enst1 = MagicMock(utr5_exons=[None] * 4, id_=1)
        enst2 = MagicMock(utr5_exons=[None] * 4, id_=2)
        enst3 = MagicMock(utr5_exons=[None] * 4, id_=3)
        enst4 = MagicMock(utr5_exons=[None] * 4, id_=4)
        enst5 = MagicMock(utr5_exons=[None] * 4, id_=5)

        selected, x, y, z = utr._utr_selection_utr5(cart, [enst1, enst2, enst3, enst4, enst5], None)
        self.assertEquals(len(selected), 1)
        self.assertEquals(selected[0].id_, 5)
        self.assertEquals(x, 'single_encompassing_utr5')
        self.assertEquals(y, '.')
        self.assertEquals(z, '.')


    @patch('cart2enst_.phase1_utr.utr_selection')
    @patch('cart2enst_.utr._enst_utr5_encompass_cart_utr5')
    def test_utr_selection_utr5_phase1_utr_selection(self, mocked_enst_utr5_encompass_cart_utr5, mocked_phase1_utr_selection):

        from selectensts_ import utr
        mocked_enst_utr5_encompass_cart_utr5.side_effect = lambda cart, enst: True
        mocked_phase1_utr_selection.side_effect = lambda candidates, log: ([candidates[1]], 'difference_type', 'decisive_criteria')

        cart = MagicMock(utr5_exons=[None] * 4)
        enst1 = MagicMock(utr5_exons=[None] * 4, id_=1)
        enst2 = MagicMock(utr5_exons=[None] * 4, id_=2)
        enst3 = MagicMock(utr5_exons=[None] * 4, id_=3)
        enst4 = MagicMock(utr5_exons=[None] * 4, id_=4)
        enst5 = MagicMock(utr5_exons=[None] * 4, id_=5)

        selected, x, y, z = utr._utr_selection_utr5(cart, [enst1, enst2, enst3, enst4, enst5], None)
        self.assertEquals(len(selected), 1)
        self.assertEquals(selected[0].id_, 2)
        self.assertEquals(x, 'phase1_utr_selection')
        self.assertEquals(y, 'difference_type')
        self.assertEquals(z, 'decisive_criteria')


    @patch('cart2enst_.utr._enst_utr5_encompass_cart_utr5')
    def test_utr_selection_utr5_no_utr5_exon_number_match(self, mocked_enst_utr5_encompass_cart_utr5):
        from selectensts_ import utr
        mocked_enst_utr5_encompass_cart_utr5.side_effect = lambda cart, enst: enst.id_ == 5

        cart = MagicMock(utr5_exons=[None] * 4)
        enst1 = MagicMock(utr5_exons=[None] * 2, id_=1)
        enst2 = MagicMock(utr5_exons=[None] * 1, id_=2)
        enst3 = MagicMock(utr5_exons=[None] * 3, id_=3)
        enst4 = MagicMock(utr5_exons=[None] * 1, id_=4)
        enst5 = MagicMock(utr5_exons=[None] * 2, id_=5)

        selected, x, y, z = utr._utr_selection_utr5(cart, [enst1, enst2, enst3, enst4, enst5], None)
        self.assertEquals(len(selected), 1)
        self.assertEquals(selected[0].id_, 5)
        self.assertEquals(x, 'single_encompassing_utr5')
        self.assertEquals(y, '.')
        self.assertEquals(z, '.')


    @patch('cart2enst_.phase1_utr.utr_selection')
    @patch('cart2enst_.utr._enst_utr5_encompass_cart_utr5')
    def test_utr_selection_no_encompassing_utr5(self, mocked_enst_utr5_encompass_cart_utr5, mocked_phase1_utr_selection):
        from selectensts_ import utr
        mocked_enst_utr5_encompass_cart_utr5.side_effect = lambda cart, enst: False
        mocked_phase1_utr_selection.side_effect =  lambda candidates, log: ([candidates[1]], 'difference_type', 'decisive_criteria')

        cart = MagicMock(utr5_exons=[None] * 4)
        enst1 = MagicMock(utr5_exons=[None] * 1, id_=1)
        enst2 = MagicMock(utr5_exons=[None] * 2, id_=2)
        enst3 = MagicMock(utr5_exons=[None] * 4, id_=3)
        enst4 = MagicMock(utr5_exons=[None] * 4, id_=4)
        enst5 = MagicMock(utr5_exons=[None] * 2, id_=5)

        selected, x, y, z = utr._utr_selection_utr5(cart, [enst1, enst2, enst3, enst4, enst5], None)
        self.assertEquals(len(selected), 1)
        self.assertEquals(selected[0].id_, 4)
        self.assertEquals(x, 'phase1_utr_selection')
        self.assertEquals(y, 'difference_type')
        self.assertEquals(z, 'decisive_criteria')


    @patch('cart2enst_.phase1_utr.utr_selection')
    @patch('cart2enst_.utr._enst_utr5_encompass_cart_utr5')
    def test_utr_selection_no_utr5_exon_number_match_and_no_encompassing_utr5(self, mocked_enst_utr5_encompass_cart_utr5, mocked_phase1_utr_selection):

        from selectensts_ import utr
        mocked_enst_utr5_encompass_cart_utr5.side_effect = lambda cart, enst: False
        mocked_phase1_utr_selection.side_effect = lambda candidates, log: ([candidates[1]], 'difference_type', 'decisive_criteria')

        cart = MagicMock(utr5_exons=[None] * 4)
        enst1 = MagicMock(utr5_exons=[None] * 1, id_=1)
        enst2 = MagicMock(utr5_exons=[None] * 2, id_=2)
        enst3 = MagicMock(utr5_exons=[None] * 3, id_=3)
        enst4 = MagicMock(utr5_exons=[None] * 5, id_=4)
        enst5 = MagicMock(utr5_exons=[None] * 2, id_=5)

        selected, x, y, z = utr._utr_selection_utr5(cart, [enst1, enst2, enst3, enst4, enst5], None)
        self.assertEquals(len(selected), 1)
        self.assertEquals(selected[0].id_, 2)
        self.assertEquals(x, 'phase1_utr_selection')
        self.assertEquals(y, 'difference_type')
        self.assertEquals(z, 'decisive_criteria')


    def test_enst_utr5_encompass_cart_utr5_forward_strand(self):

        from selectensts_ import utr

        transcript_1 = MagicMock(transcript_start=1000, transcript_end=2000, strand='+')
        transcript_2 = MagicMock(transcript_start=900, transcript_end=1700, strand='+')
        transcript_3 = MagicMock(transcript_start=1000, transcript_end=2500, strand='+')
        transcript_4 = MagicMock(transcript_start=1200, transcript_end=2000, strand='+')

        self.assertTrue(utr._enst_utr5_encompass_cart_utr5(transcript_1, transcript_2))
        self.assertTrue(utr._enst_utr5_encompass_cart_utr5(transcript_1, transcript_3))
        self.assertFalse(utr._enst_utr5_encompass_cart_utr5(transcript_1, transcript_4))
        self.assertTrue(utr._enst_utr5_encompass_cart_utr5(transcript_4, transcript_4))


    def test_enst_utr5_encompass_cart_utr5_reverse_strand(self):

        from selectensts_ import utr

        transcript_1 = MagicMock(transcript_start=1000, transcript_end=2000, strand='-')
        transcript_2 = MagicMock(transcript_start=1200, transcript_end=1700, strand='-')
        transcript_3 = MagicMock(transcript_start=1100, transcript_end=2500, strand='-')
        transcript_4 = MagicMock(transcript_start=900, transcript_end=2000, strand='-')

        self.assertFalse(utr._enst_utr5_encompass_cart_utr5(transcript_1, transcript_2))
        self.assertTrue(utr._enst_utr5_encompass_cart_utr5(transcript_1, transcript_3))
        self.assertTrue(utr._enst_utr5_encompass_cart_utr5(transcript_1, transcript_4))
        self.assertTrue(utr._enst_utr5_encompass_cart_utr5(transcript_4, transcript_4))


class TestUTRDifferences(TestCase):


    def setUp(self):

        self.transcript_1 = MagicMock(utr5_exons=[(100, 200), (300, 400)], utr3_exons=[(1000, 1200)])
        self.transcript_2 = MagicMock(utr5_exons=[(300, 400)], utr3_exons=[(1000, 1200)])
        self.transcript_3 = MagicMock(utr5_exons=[(100, 200), (300, 400)], utr3_exons=[(1020, 1200)])
        self.transcript_4 = MagicMock(utr5_exons=[(100, 200), (300, 400)], utr3_exons=[])
        self.transcript_5 = MagicMock(utr5_exons=[], utr3_exons=[])


    def test_utr_difference(self):

        from selectensts_ import utr

        self.assertEquals(utr.utr_difference(self.transcript_1, self.transcript_1), '.')
        self.assertEquals(utr.utr_difference(self.transcript_1, self.transcript_2), 'UTR5')
        self.assertEquals(utr.utr_difference(self.transcript_1, self.transcript_3), 'UTR3')
        self.assertEquals(utr.utr_difference(self.transcript_2, self.transcript_3), 'UTR5,UTR3')
        self.assertEquals(utr.utr_difference(self.transcript_1, self.transcript_4), 'UTR3')
        self.assertEquals(utr.utr_difference(self.transcript_4, self.transcript_5), 'UTR5')


    def test_utr_exon_number_difference(self):

        from selectensts_ import utr

        self.assertEquals(utr.utr_exon_number_difference(self.transcript_1, self.transcript_2), 'UTR5:-1')
        self.assertEquals(utr.utr_exon_number_difference(self.transcript_2, self.transcript_1), 'UTR5:+1')
        self.assertEquals(utr.utr_exon_number_difference(self.transcript_4, self.transcript_5), 'UTR5:-2')
        self.assertEquals(utr.utr_exon_number_difference(self.transcript_3, self.transcript_5), 'UTR5:-2,UTR3:-1')
