#pragma once
#ifndef IDOCK_KERNEL_HPP
#define IDOCK_KERNEL_HPP

class kernel
{
protected:
	kernel(const int device_id) : device_id(device_id)
	{
	}

	const int device_id;
};

#endif
