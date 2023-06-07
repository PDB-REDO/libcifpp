#include <cif++.hpp>

class dummy_parser : public cif::sac_parser
{
  public:
	dummy_parser(std::istream &is)
		: sac_parser(is)
	{
	}

	void produce_datablock(std::string_view name) override
	{
	}

	void produce_category(std::string_view name) override
	{
	}

	void produce_row() override
	{
	}

	void produce_item(std::string_view category, std::string_view item, std::string_view value) override
	{
	}
};


int main()
{
	cif::gzio::ifstream in("/srv/data/pdb/mmCIF/gl/8glv.cif.gz");

	dummy_parser parser(in);
	parser.parse_file();

	// cif::file f("/srv/data/pdb/mmCIF/gl/8glv.cif.gz");

	return 0;
}