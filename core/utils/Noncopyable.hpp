#ifndef DMP_UTILS_NONCOPYABLE_HPP
#define DMP_UTILS_NONCOPYABLE_HPP

namespace dmp
{
	class Noncopyable
	{
	public:
		Noncopyable( Noncopyable const & ) = delete;
		Noncopyable& operator=( Noncopyable const & ) = delete;

		Noncopyable() = default;
	};
}

#endif //DMP_UTILS_NONCOPYABLE_HPP
