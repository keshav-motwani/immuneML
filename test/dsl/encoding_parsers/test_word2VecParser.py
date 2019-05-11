from unittest import TestCase

from source.dsl.encoding_parsers.Word2VecParser import Word2VecParser
from source.encodings.word2vec.model_creator.ModelType import ModelType


class TestWord2VecParser(TestCase):
    def test_parse(self):
        model = {
            "k": 3,
            "size": 16
        }

        parsed, specs = Word2VecParser.parse(model)
        self.assertEqual(3, parsed["k"])
        self.assertEqual(ModelType.SEQUENCE, parsed["model_creator"])
        self.assertEqual(16, parsed["size"])

        parsed, specs = Word2VecParser.parse({})
        self.assertEqual(3, parsed["k"])
        self.assertEqual(ModelType.SEQUENCE, parsed["model_creator"])
        self.assertEqual(64, parsed["size"])

        model["k"] = 0
        self.assertRaises(AssertionError, Word2VecParser.parse, model)