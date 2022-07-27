	// void write(std::ostream &os, const std::vector<size_t> &order, bool includeEmptyColumns)
	// {
	// 	if (empty())
	// 		return;

	// 	// If there are multiple rows in this category, we need a _loop
	// 	if (size() == 1)
	// 	{
	// 		os << "loop_" << std::endl;

	// 		std::vector<size_t> columnWidths;

	// 		for (auto cix : order)
	// 		{
	// 			auto &col = mColumns[cix];
	// 			os << '_' << mName << '.' << col.mName << ' ' << std::endl;
	// 			columnWidths.push_back(2);
	// 		}

	// 		for (auto Row = mHead; Row != nullptr; Row = Row->mNext)
	// 		{
	// 			for (auto v = Row->mValues; v != nullptr; v = v->mNext)
	// 			{
	// 				if (strchr(v->mText, '\n') == nullptr)
	// 				{
	// 					size_t l = strlen(v->mText);

	// 					if (not isUnquotedString(v->mText))
	// 						l += 2;

	// 					if (l > 132)
	// 						continue;

	// 					if (columnWidths[v->mColumnIndex] < l + 1)
	// 						columnWidths[v->mColumnIndex] = l + 1;
	// 				}
	// 			}
	// 		}

	// 		for (auto Row = mHead; Row != nullptr; Row = Row->mNext) // loop over rows
	// 		{
	// 			size_t offset = 0;

	// 			for (size_t cix : order)
	// 			{
	// 				size_t w = columnWidths[cix];

	// 				std::string s;
	// 				for (auto iv = Row->mValues; iv != nullptr; iv = iv->mNext)
	// 				{
	// 					if (iv->mColumnIndex == cix)
	// 					{
	// 						s = iv->mText;
	// 						break;
	// 					}
	// 				}

	// 				if (s.empty())
	// 					s = "?";

	// 				size_t l = s.length();
	// 				if (not isUnquotedString(s.c_str()))
	// 					l += 2;
	// 				if (l < w)
	// 					l = w;

	// 				if (offset + l > 132 and offset > 0)
	// 				{
	// 					os << std::endl;
	// 					offset = 0;
	// 				}

	// 				offset = detail::writeValue(os, s, offset, w);

	// 				if (offset > 132)
	// 				{
	// 					os << std::endl;
	// 					offset = 0;
	// 				}
	// 			}

	// 			if (offset > 0)
	// 				os << std::endl;
	// 		}
	// 	}
	// 	else
	// 	{
	// 		// first find the indent level
	// 		size_t l = 0;

	// 		for (auto &col : mColumns)
	// 		{
	// 			std::string tag = '_' + mName + '.' + col.mName;

	// 			if (l < tag.length())
	// 				l = tag.length();
	// 		}

	// 		l += 3;

	// 		for (size_t cix : order)
	// 		{
	// 			auto &col = mColumns[cix];

	// 			os << '_' << mName << '.' << col.mName << std::string(l - col.mName.length() - mName.length() - 2, ' ');

	// 			std::string s;
	// 			for (auto iv = mHead->mValues; iv != nullptr; iv = iv->mNext)
	// 			{
	// 				if (iv->mColumnIndex == cix)
	// 				{
	// 					s = iv->mText;
	// 					break;
	// 				}
	// 			}

	// 			if (s.empty())
	// 				s = "?";

	// 			size_t offset = l;
	// 			if (s.length() + l >= kMaxLineLength)
	// 			{
	// 				os << std::endl;
	// 				offset = 0;
	// 			}

	// 			if (detail::writeValue(os, s, offset, 1) != 0)
	// 				os << std::endl;
	// 		}
	// 	}

	// 	os << "# " << std::endl;
	// }

	void write(std::ostream &os) const
	{
		// std::vector<size_t> order(mColumns.size());
		// iota(order.begin(), order.end(), 0);
		// write(os, order, false);

		os << '#' << m_name << std::endl;
		for (auto &r : *this)
		{
			for (auto &f : r)
				os << '_' << m_name << '.' << f.name() << ' ' << f.value() << std::endl;
		}
	}

	// void Category::write(std::ostream &os, const std::vector<std::string> &columns)
	// {
	// 	// make sure all columns are present
	// 	for (auto &c : columns)
	// 		addColumn(c);

	// 	std::vector<size_t> order;
	// 	order.reserve(mColumns.size());

	// 	for (auto &c : columns)
	// 		order.push_back(getColumnIndex(c));

	// 	for (size_t i = 0; i < mColumns.size(); ++i)
	// 	{
	// 		if (std::find(order.begin(), order.end(), i) == order.end())
	// 			order.push_back(i);
	// 	}

	// 	write(os, order, true);
	// }
